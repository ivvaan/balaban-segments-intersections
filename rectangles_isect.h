#pragma once
/*
*
*      Copyright (c)  2011-2026  Ivan Balaban
*      ivvaan@gmail.com
*
This file is part of Seg_int library.

Seg_int is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Seg_int is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Seg_int.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "collection_base.h"
#include "utils.h"  
#include <iostream>
#include <cassert>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <bit>

namespace SegmentTreeAndList {
  // 3 or 4 is the best value for P, 
  // 3 is better for small n, because it has smaller constant factors, but 4 is better for large n, 
  // because it has smaller memory usage and more cache friendly
  // so we use 4 by default, but you can change it to 3 if you want.
  template<int4 P = 4> //2^P -ary tree is used to store filling info and search for previous filled element in list, P must be at least 1
  struct TreeList {
    constexpr static const int4 S = (1 << P) - 1;
    constexpr static const int4 Q = S - 1;
    constexpr static const int4 M = ~S;
    int4* tree_list = nullptr;//[0..sz-1] - used by binary tree to store filling subtree info, 
    //[sz..SZ-1] - list of filled elements, tree_list[i] is the next filled element in list after i, 
    // tree_list[0] is the head of list (the first filled element in list or 1 if list is empty)
    // if list is not empty, tree_list[i] is always >= 1 
    int4 n = 0;
    int4 sz = 0;//first power of 2 greater or equal to n
    int4 SZ = 0;//sz+n

    bool is_filled(int4 pos) const {
      return tree_list[pos];
    }

    static int4 get_rightmost_son(int4 father) {
      return 1 + (father << P);
    }

    static int4 get_leftmost_sibling(int4 pos) {
      return ((pos + Q) & M) - Q;
    }

    int4 get_left_filled_sibling(int4 pos) const {
      auto ls = get_leftmost_sibling(pos);
      --pos;
      while (ls <= pos) {
        if (is_filled(pos))
          return pos;
        --pos;
      }
      return 0;
    };

    int4 get_rightmost_filled_son(int4 father) const {
      int4 i = get_rightmost_son(father);
      while (!is_filled(i))
        --i;
      return i;
    };

    static int4 get_father(int4 pos) {
      return (pos + Q) >> P;
    }

    static constexpr const bool lst_del = false;
    static constexpr const bool lst_ins = true;

    template<bool is_insert>
    void list_change(int4 rank) {
      auto pos = sz + rank;
      int4 prev_elem = 0;
      while (pos != 1) {
        if (prev_elem = get_left_filled_sibling(pos)) {
          pos = get_father(pos);
          break;
        }
        if constexpr (is_insert) {
          ++tree_list[pos = get_father(pos)];
        }
        else {
          --tree_list[pos = get_father(pos)];
          assert(tree_list[pos] >= 0);
        }
      };

      if (prev_elem) {
        while (pos) {
          if constexpr (is_insert) {
            ++tree_list[pos];
          }
          else {
            --tree_list[pos];
            assert(tree_list[pos] >= 0);
          }
          pos = get_father(pos);
        }
        while (prev_elem < sz)
          prev_elem = get_rightmost_filled_son(prev_elem);
      }
      pos = sz + rank;

      if constexpr (is_insert) {
        tree_list[pos] = tree_list[prev_elem];
        tree_list[prev_elem] = pos;
      }
      else {
        assert(tree_list[prev_elem] == pos);
        tree_list[prev_elem] = tree_list[pos];
        tree_list[pos] = 0;
      }
    }

    void list_insert(int4 rank) {
      list_change<lst_ins>(rank);
    }

    void list_delete(int4 rank) {
      list_change<lst_del>(rank);
    }

    int4 get_next(int4 rank) const {
      return tree_list[sz + rank] - sz;
    }

    int4 get_sz() {
      int4 v = 1;
      int4 r = 1;
      while (v < n) {
        r += v;
        v <<= P;
      }
      return r;
    }

    TreeList(int4 _n) : n(_n) {
      sz = get_sz();
      SZ = sz + n;
      tree_list = new int4[SZ];
      std::fill_n(tree_list, SZ, 0);
      tree_list[0] = 1;//fake last list element is stored in 0 position - the header of list
      //ordered_list[i]==1 where i>=sz means that i is filled, but next filled element is not exists (i.e.last element in list is i)
    }

    ~TreeList() {
      if (tree_list) {
        delete[] tree_list;
        tree_list = nullptr;
      }
    }

  };

  template<>
  int4 TreeList<1>::get_left_filled_sibling(int4 pos) const {
    if ((pos & 1) == 0) return 0;
    auto sibling = pos ^ 1;
    return is_filled(sibling) ? sibling : 0;
  }

  template<>
  int4 TreeList<1>::get_rightmost_filled_son(int4 pos) const {
    auto rs = get_rightmost_son(pos);
    return is_filled(rs) ? rs : rs - 1;
  }


  struct STree {
    int4 n = 0;
    int4 sz = 0;//first power of 2 greater or equal to n
    int4 SZ = 0;//sz+n
    int4 depth = 0;

    int4 get_tree_size() const {
      return SZ >> 1;
    }

    template<typename action_func>
    void locate(int4 rank, int4 id, action_func do_locate) {
      auto pos = (sz + rank) >> 1;
      do {
        do_locate(pos, id);
        pos >>= 1;
      } while (pos > 0);
    }

    template<typename action_func>
    void insert_range(int4 l, int4 r, int4 id, action_func do_insert) {
      if ((l == 0) && (r == n)) {
        do_insert(1, id);
        return;
      }
      int4 pow = depth;
      int4 from, to;
      do {
        --pow;
        from = (l >> pow);
        from += l != (from << pow);
        to = r >> pow;
      } while (from >= to);
      if (pow < 1) return;
      auto idx = (sz >> pow) + from;
      do_insert(idx, id);
      if (from + 1 < to)
        do_insert(idx + 1, id);
      if (pow < 2) return;
      insert_left(l, from << pow, id, do_insert);
      insert_right(to << pow, r, id, do_insert);
    }

    template<typename action_func>
    void insert_left(int4 l, int4 r, int4 id, action_func do_insert) {
      if (l + 2 > r)
        return;
      auto l_next = l + 1;
      do {
        auto pow = std::bit_width((unsigned)(r - l)) - 1;
        r -= (1 << pow);
        do_insert((sz + r) >> pow, id);
      } while (l_next < r);
    }

    template<typename action_func>
    void insert_right(int4 l, int4 r, int4 id, action_func do_insert) {
      if (l + 2 > r)
        return;
      auto r_prev = r - 1;
      do {
        auto pow = std::bit_width((unsigned)(r - l)) - 1;
        do_insert((sz + l) >> pow, id);
        l += (1 << pow);
      } while (l < r_prev);
    }

    int4 get_sz() {
      int4 v = 1;
      while (v < n) v <<= 1;
      return v;
    }

    int4 get_depth() {
      int4 d = 0;
      for (int4 v = n; v > 0; v >>= 1)
        ++d;
      return d;
    }

    STree(int4 _n) : n(_n) {
      sz = get_sz();
      SZ = sz + n;
      depth = get_depth();
    }

    ~STree() {}

  };

}

template<typename SegmCollection>
void rect_find_intersections(SegmCollection &coll) {

  struct count_info {
    int4 ins = 0, del = 0, max = 0;
    void insert() {if (++ins > max) max = ins;}
    void erase() {++del;}
    void locate() {ins -= del;del = 0;}
  };

  struct arr_info {
    int4 beg = 0, end = 0;
  };

  using namespace SegmentTreeAndList;
  int4 n = coll.GetSegmNumb();

  auto ranks2pointsX_keeper = std::make_unique<int[]>(n * 2);
  auto ranks2pointsX = coll.get_sorted_bounds(ranks2pointsX_keeper.get(), false);
  auto ranks2pointsY_keeper = std::make_unique<int[]>(n * 2);
  auto ranks2pointsY = coll.get_sorted_bounds(ranks2pointsY_keeper.get(), true);

  auto points2ranksY_keeper = std::make_unique<int[]>(n * 2);
  auto points2ranksY = points2ranksY_keeper.get();
  for (int i = 0; i < n * 2; ++i) {
    points2ranksY[ranks2pointsY[i]] = i;
  }
  auto ranks2rectanglesY = ranks2pointsY;//we dont need ranks2pointsY, but we need ranks2rectanglesY, so we can reuse allocated array
  for (int i = 0; i < n * 2; ++i) {
    ranks2rectanglesY[i] >>= 1;//rectangle id is any endpoint id divided by 2, because each rectangle has 2 endpoints
  }

  TreeList<> tl(n * 2);
  STree st(n * 2);
  //first we need to calculate how many rectangles can be in each node of segment tree at the same time, to allocate arrays for nodes
  auto tree_size = st.get_tree_size();
  auto counts_keeper = std::make_unique<count_info[]>(tree_size);
  auto counts = counts_keeper.get();
  std::fill_n(counts, tree_size, count_info{});

  for (int i = 0; i < n * 2; ++i) {
    auto point = ranks2pointsX[i];
    auto rect_id = point >> 1;
    if ((point & 1) == 0) {
      st.locate(points2ranksY[rect_id << 1], rect_id, [=](int pos, int id) {
        assert(pos < tree_size);
        counts[pos].locate();
        });
      st.insert_range(points2ranksY[rect_id << 1], points2ranksY[(rect_id << 1) + 1], rect_id, [=](int pos, int id) {
        assert(pos < tree_size);
        counts[pos].insert();
        });
    }
    else {
      st.insert_range(points2ranksY[rect_id << 1], points2ranksY[(rect_id << 1) + 1], rect_id, [=](int pos, int id) {
        assert(pos < tree_size);
        counts[pos].erase();
        });
    }
  }
  //now we know the maximum count of rectangles in each node, we can allocate arrays for nodes
  auto node_arrays_keeper = std::make_unique<arr_info[]>(tree_size);
  auto node_arrays = node_arrays_keeper.get();
  int acc = 0;
  for (int i = 0; i < tree_size; ++i) {
    node_arrays[i].beg = node_arrays[i].end = acc;
    acc += counts[i].max;
  }
  counts_keeper.reset();
  // With per-node arrays allocated, run the main algorithm: sweep along the X axis.
  // For each rectangle, locate its lower Y endpoint in the segment tree to find
  // all rectangles that overlap it in Y. Then insert the rectangle for future checks
  // and scan the corner list between its Y endpoints, since those also intersect.
  // Duplicates can appear both in tree nodes and the corner list, so use
  // `dublicate_checker` to avoid reporting the same pair more than once.
  auto node_rects_keeper = std::make_unique<int4[]>(acc);
  auto node_rects = node_rects_keeper.get();
  auto is_removed_keeper = std::make_unique<char[]>(n);
  auto is_removed = is_removed_keeper.get();
  std::fill_n(is_removed, n, 0);
  auto dublicate_checker_keeper = std::make_unique<int4[]>(n);
  auto dublicate_checker = dublicate_checker_keeper.get();
  std::fill_n(dublicate_checker, n, -1);
  for (int4 i = 0; i < n * 2; ++i) {
    auto point = ranks2pointsX[i];
    auto rect_id = point >> 1;
    coll.SetCurSeg(rect_id);
    auto beginY = points2ranksY[rect_id << 1];
    auto endY = points2ranksY[point | 1];

    if ((point & 1) == 0) {
      auto do_locate = [=,&coll](int4 pos, int4 rect) {
        auto& arr_info = node_arrays[pos];
        int4 new_end = arr_info.beg;
        for (int4 i = arr_info.beg; i < arr_info.end; ++i) {
          auto other = node_rects[i];
          //lazy deletion - we mark rectangles as removed, but do not remove them from arrays, 
          // until we need to locate new rectangle, then we do remove all removed rectangles 
          // from array and update end index.
          if (is_removed[other])
            continue;
          node_rects[new_end++] = other;
          dublicate_checker[other] = rect;
          coll.SSCurSegIntWith(other);
        }
        arr_info.end = new_end;
        };
      st.locate(beginY, rect_id, do_locate);

      auto do_insert = [=](int4 pos, int4 id) {
        node_rects[node_arrays[pos].end++] = id;
        };
      st.insert_range(beginY, endY, rect_id, do_insert);
      tl.list_insert(endY);
      tl.list_insert(beginY);
      auto next = tl.get_next(beginY);
      while (next != endY) {
        auto other_id = ranks2rectanglesY[next];
        next = tl.get_next(next);
        if (dublicate_checker[other_id] == rect_id)
          continue;
        dublicate_checker[other_id] = rect_id;
        coll.SSCurSegIntWith(other_id);
      }
    }
    else {
      tl.list_delete(beginY);
      tl.list_delete(endY);
      //we do not remove rectangle from arrays, but mark it as removed, so it will be removed later when we need to locate new rectangle
      is_removed[rect_id] = 1;
    }
  }
}
