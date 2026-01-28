

class CRemaper {
  struct MarkedListBounds {
    uint4 beg = 0;
    uint4 end = 0;
    uint4 flags = 0;
    bool is_mapped_entirely() const {
      return (flags != 0);
    };
  };

  static REAL get_rot_angle(uint4 n, TLineSegment1 sc[]) {
    REAL* arr;
    DECL_RAII_ARR(arr, n + 2);
    arr[0] = -M_PI / 2.;
    std::transform(sc, sc + n, arr,
      [](const TLineSegment1& s) {
        return std::atan2(s.shift.y, s.shift.x);
      });
    arr[n + 1] = M_PI / 2.;
    auto max_gap_mid = get_max_gap_middle(n + 2, arr);
    return max_gap_mid + (max_gap_mid > 0 ? -M_PI / 2. : M_PI / 2.);
  }

public:

  template<bool remove_zero_seg>
  uint4 int_seg_from_real(uint4 n, TLineSegment1 sc[], std::vector<TIntegerVect>& points) {
    //REAL angle = CRandomValueGen().GetRandomDouble();
    //REAL angle = get_rot_angle(n,sc);//try to rotate the scene to remove vertical segments
    //REAL si = std::sin(angle), co = std::cos(angle);

    REAL si = 0., co = 1.;
    auto mm_rect = get_rot_minmax(n, sc, si, co);
    auto x_transf = transf1D{ .shift = -mm_rect.ld.x - mm_rect.get_width() / 2.,
      .scale = 2. * range / mm_rect.get_width() };
    auto y_transf = transf1D{ .shift = -mm_rect.ld.y - mm_rect.get_height() / 2.,
      .scale = 2. * range / mm_rect.get_height() };
    seg = set_size(seg_v, n);
#ifdef PRINT_SEG_AND_INT
    TIntegerSegment::coll_begin = coll_begin;
#endif
    TIntegerSegment s;
    uint4 i = 0;
    for (uint4 m = 0; m != n; ++m) {
      s.Init(sc[m], si, co, x_transf, y_transf);
      if constexpr (remove_zero_seg) {
        if (s.shift.is_zero()) continue;
      }
      seg[i] = s;
      points[first_point(i)] = s.BegPoint();
      points[last_point(i)] = s.EndPoint();
      ++i;
    }
    return i;
  }

  bool prepare_remap(std::vector<uint4>& indexes, TIntegerVect points[],
    std::vector<TIntegerSegment>& res_coll,
    std::vector<TIntegerVect>& res_pts) {

    auto N = indexes.size();
    uint4 remap_size = 0;
    int seg_numb = 0;
    remapped_SN = 0;
    for (auto pt : indexes) {
      seg_numb += is_last(pt) ? -1 : 1;
      assert(seg_numb >= 0);
      remap_size += seg_numb;
      remapped_SN += (seg_numb) ? 1 : 0;
    }
    if ((initial_SN == remapped_SN) && (initial_SN == nonzero_N))
      return true;//not remapped
    rrec = set_size(rrec_v, remapped_SN);
    auto remap = set_size(remap_v, remap_size);

    set_size(res_coll, remapped_SN);
    set_size(res_pts, 2 * remapped_SN);
    uint4 remaper_pos = 0, new_seg_num = 0;
    std::vector<uint4> stack;
    stack.reserve(nonzero_N);
    bool not_remapped = initial_SN == nonzero_N;
    for (size_t i = 0; i < N; ++i) {
      auto pt = indexes[i];
      if (is_last(pt)) {
        remove_by_val(stack, get_segm(pt));
      }
      else {
        auto s = get_segm(pt);
        for (auto p : stack) {
          registrator->register_pair(s, p);
        }
        stack.push_back(s);
      };
      auto size = stack.size();

      if (size != 0) {
        assert(i + 1 < N);

        auto next_pt = indexes[i + 1];
        auto beg = points[pt];
        auto end = points[next_pt];

        if (beg != end) {
          std::copy(stack.begin(), stack.end(), remap + remaper_pos);
          uint4 is_mapped_entirely = 0;
          if (size == 1) {
            auto s = remap[remaper_pos];
            if (pt == first_point(s) && (next_pt == last_point(s)))
              is_mapped_entirely = 1;
          }
          rrec[new_seg_num] = { remaper_pos,remaper_pos + (uint4)size, is_mapped_entirely };
          remaper_pos += size;
          not_remapped &= size < 2;
          res_coll[new_seg_num] = { beg,end };
          res_pts[first_point(new_seg_num)] = beg;
          res_pts[last_point(new_seg_num)] = end;
          ++new_seg_num;

        }
      };
    }
    remapped_SN = new_seg_num;
    remap_size = remaper_pos;
    return not_remapped;
  }

  void register_pair(uint4 s1, uint4 s2, uint4 int_type) {
    auto orig_int_type = seg[s1].get_int_type_beg(seg[s2]);
    if (int_type == orig_int_type)
      registrator->register_pair(s1, s2);
  }

  void register_pair(uint4 s1, uint4 s2) {
    auto rr1 = rrec[s1];
    auto rr2 = rrec[s2];

    if (rr1.is_mapped_entirely() && rr2.is_mapped_entirely()) {
      //if we have full mapping segment to segment, no need to check anything
      //trivial remap_v (should be most of the time)
      registrator->register_pair(remap_v[rr1.beg], remap_v[rr2.beg]);
      return;
    }

    auto int_type = remapped_segs[s1].get_int_type_beg(remapped_segs[s2]);

    if (int_type == _IntType::common_int) {
      //if nondegenerate common intersection of mapped segments s1 and s2 -> 
      //trivial remap_v for both one to one or many to many mapping

      if (rr1.end + rr2.end - rr1.beg - rr2.beg < 3) {//one to one mapping 
        registrator->register_pair(remap_v[rr1.beg], remap_v[rr2.beg]);
        return;
      };

      for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1) {//many to many mapping
        for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
          registrator->register_pair(remap_v[rr1.beg], remap_v[rr2.beg]);
      };
      return;
    };
    //non common intersection, need to check more carefully
    //need to check type of intersection to exclude double registration

    if (rr1.end + rr2.end - rr1.beg - rr2.beg < 3) {
      register_pair(remap_v[rr1.beg], remap_v[rr2.beg], int_type);
      return;
    };

    //most complex case, need to check for nonunique segments; hopefully rare
    // 
    //marking all nonunique segments
    constexpr const uint4 mark = 1 << (sizeof(uint4) * 8 - 1);
    for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
      for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
        if (remap_v[seg1] == remap_v[seg2]) {//same segment found, mark it 
          remap_v[seg1] |= mark;
          remap_v[seg2] |= mark;
          break;
        }

    //register only intersections with nonmarked segments, otherwise we can have double registration
    for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
      if ((remap_v[seg1] & mark) == 0) {
        for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
          if ((remap_v[seg2] & mark) == 0)
            //also need to check type of intersection to exclude double registration
            register_pair(remap_v[seg1], remap_v[seg2], int_type);
      }

    // erase marks, back to original values
    for (uint4 seg1 = rr1.beg; seg1 < rr1.end; ++seg1)
      remap_v[seg1] &= ~mark;
    for (uint4 seg2 = rr2.beg; seg2 < rr2.end; ++seg2)
      remap_v[seg2] &= ~mark;

  };

  template<bool remove_zero_seg = true>
  auto NoRemapInit(uint4 n, TLineSegment1* sc, IntersectionRegistrator* r, int4 range, CTHIS& collection) {
    initial_SN = n;
    std::vector<TIntegerVect> points(2 * n);

    this->range = range;
    remapped_SN = nonzero_N = int_seg_from_real<remove_zero_seg>(n, sc, points);
    if constexpr (remove_zero_seg)
      initial_SN = nonzero_N;
    registrator = r;

    collection.segments = std::move(seg_v);
    collection.points = std::move(points);
    return;
  };

  template<bool remove_zero_seg = true>
  auto NoRemapInit(uint4 n, TIntegerSegment* sc, IntersectionRegistrator* r, int4 range, CTHIS& collection) {
    initial_SN = n;
    std::vector<TIntegerVect> points(2 * n);

    this->range = range;
    uint4 i = 0;
    seg = set_size(seg_v, n);
    for (uint4 m = 0; m != n; ++m) {
      if constexpr (remove_zero_seg) {
        if (sc[m].shift.is_zero()) continue;
      }
      seg[i] = sc[m];
      points[first_point(i)] = sc[m].BegPoint();
      points[last_point(i)] = sc[m].EndPoint();
      ++i;
    }
    remapped_SN = nonzero_N = i;
    if constexpr (remove_zero_seg)
      initial_SN = nonzero_N;
    registrator = r;

    collection.segments = std::move(seg_v);
    collection.points = std::move(points);
    return;
  };

  template<bool remove_zero_seg = true>
  auto FromIntSegVect(std::vector<TIntegerSegment>& v, IntersectionRegistrator* r, CTHIS& collection) {
    auto n = nonzero_N = v.size();
    std::vector<TIntegerVect> points(2 * n);
    uint4 i = 0;
    for (auto& s : v) {
      if (s.shift.is_zero()) {
        --nonzero_N;
        continue;
      }
      assert(s.is_inside_int_range());
      seg_v.push_back(s);
      points[first_point(i)] = s.BegPoint();
      points[last_point(i)] = s.EndPoint();
      ++i;
    }
    seg = seg_v.data();
    remapped_SN = nonzero_N;
    if constexpr (remove_zero_seg)
      initial_SN = nonzero_N;
    registrator = r;

    collection.segments = std::move(seg_v);
    collection.points = std::move(points);
    return;
  };

  auto TurnRemapOn(CTHIS& collection) {
    auto points = std::move(collection.points);
    seg_v = std::move(collection.segments);
    seg = seg_v.data();
    std::vector<uint4> indexes;
    auto n = initial_SN;
    indexes.reserve(2 * n);

    for (uint4 i = 0; i < 2 * n; ++i)
      if (seg[get_segm(i)].shift.is_non_zero())
        indexes.push_back(i);

    nonzero_N = indexes.size() / 2;

    auto comparator = [segments = seg, pts = points.data()](uint4 pt1, uint4 pt2) {
      auto i1 = get_segm(pt1);
      auto i2 = get_segm(pt2);
      if (i1 == i2)
        return pt1 < pt2; //same segment
      auto& s1 = segments[i1];
      auto& s2 = segments[i2];
      auto S = s1.shift % s2.shift;
      if (S != 0)//segments non parallel
        return S < 0;
      //segments parallel
      auto shift = s1.shift;// +s2.shift;
      assert(shift.is_non_zero());
      auto oo = pts[pt2] - pts[pt1];
      S = shift % oo;
      if (S != 0)// segments on different lines
        return  S < 0;
      //segments parallel and lies on one line
      auto prod = oo * shift;
      if (prod != 0)// points not coinside
        return  0 < prod;
      // points coinside
      if (is_last(pt1) != is_last(pt2))
        return is_last(pt1) > is_last(pt2); //different segments end first
      return pt1 < pt2;
      };

    std::sort(indexes.begin(), indexes.end(), comparator);

    bool not_remapped = prepare_remap(indexes, points.data(), collection.segments, collection.points);
    if (not_remapped) {
      collection.segments = std::move(seg_v);
      collection.points = std::move(points);
      return false;
    }
    remapped_segs = collection.segments.data();
    remapped_ends = collection.points.data();

    return true;
  };

  TIntegerSegment& get_original_segment(uint4 s) {
    return seg[s];
  }


  auto get_N() {
    return remapped_SN;
  }

  IntersectionRegistrator* registrator = nullptr;

  void clone_from(CRemaper* other) {
    initial_SN = other->initial_SN;
    nonzero_N = other->nonzero_N;
    remapped_SN = other->remapped_SN;
    remap_size = other->remap_size;
    range = other->range;
    registrator = other->registrator;
    remap_v = other->remap_v;//the data mast be copied when cloning (read\write access)
    //just copy pointers (for read only access)
    rrec = other->rrec;
    seg = other->seg;
    remapped_segs = other->remapped_segs;
    remapped_ends = other->remapped_ends;
  };

  void unclone() {
    rrec = nullptr;
    remapped_segs = nullptr;
    remapped_ends = nullptr;
    registrator = nullptr;
  };

private:
  uint4 initial_SN = 0;
  uint4 nonzero_N = 0;
  uint4 remapped_SN = 0;
  uint4 remap_size = 0;
  int4 range = 1;
  std::vector<MarkedListBounds> rrec_v;
  MarkedListBounds* rrec = nullptr;
  std::vector<uint4> remap_v;//mast be copied when cloning
  std::vector<TIntegerSegment> seg_v;//if remapped original int segments are stored here
  TIntegerSegment* seg = nullptr;
  TIntegerSegment* remapped_segs = nullptr;
  TIntegerVect* remapped_ends = nullptr;
  //CRemaper* clone_of = nullptr;

};


