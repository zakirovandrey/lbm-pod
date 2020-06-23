inline __device__ ftype Hermite(const int ni, const ftype3 _v3, const int abc[]) {
  const ftype v[3] = {_v3.x,_v3.y,_v3.z};
  const int a=abc[0], b=abc[1], c=abc[2], d=abc[3], e=abc[4];
  if(ni==0) return 1; else
  if(ni==1) return v[a]; else
  if(ni==2) return v[a]*v[b] - (a==b); else
  if(ni==3) return v[a]*v[b]*v[c] - ( v[a]*(b==c) + v[b]*(a==c) + v[c]*(a==b) ); else
  if(ni==4) return v[a]*v[b]*v[c]*v[d] - ( v[a]*v[b]*(c==d)+ v[a]*v[c]*(b==d)+ v[a]*v[d]*(b==c)+ v[b]*v[c]*(a==d)+ v[b]*v[d]*(a==c)+ v[c]*v[d]*(a==b) )+
                                     ( (a==b)*(c==d) + (a==c)*(b==d) + (a==d)*(b==c) );    else
  if(ni==5) return v[a]*v[b]*v[c]*v[d]*v[e] 
                     - ( v[a]*v[b]*v[c]*(d==e) + v[a]*v[b]*v[d]*(c==e) + v[a]*v[b]*v[e]*(c==d) + v[a]*v[c]*v[d]*(b==e) + v[a]*v[c]*v[e]*(b==d) + v[a]*v[d]*v[e]*(b==c)
                                               + v[b]*v[c]*v[d]*(a==e) + v[b]*v[c]*v[e]*(a==d) + v[b]*v[d]*v[e]*(a==c) + v[c]*v[d]*v[e]*(a==b) )
                     - ( v[a]*(b==c)*(d==e) + v[a]*(b==d)*(c==e) + v[a]*(b==e)*(c==d) +
                         v[b]*(a==c)*(d==e) + v[b]*(a==d)*(c==e) + v[b]*(a==e)*(c==d) + 
                         v[c]*(a==b)*(d==e) + v[c]*(a==d)*(b==e) + v[c]*(a==e)*(b==d) + 
                         v[d]*(a==b)*(c==e) + v[d]*(a==c)*(b==e) + v[d]*(a==e)*(b==c) + 
                         v[e]*(a==b)*(c==d) + v[e]*(a==c)*(b==d) + v[e]*(a==d)*(b==c)     ) 
            ; else
  assert(0);
  return 0;
}

template<int _ORD> inline constexpr int get_all_tensor_coeffs(){ return get_power<DIM,_ORD>::value + get_all_tensor_coeffs<_ORD-1>(); }
template<> inline constexpr int get_all_tensor_coeffs<-1>(){ return 0; }
template<int MAX_ORDER> struct TensorCoeffs{
  static const int Order=MAX_ORDER;
  static const int NCoeffs = get_all_tensor_coeffs<MAX_ORDER>();
  ftype k[NCoeffs?NCoeffs:1];
  __host__ __device__ inline TensorCoeffs() {}
  __host__ __device__ inline TensorCoeffs(const ftype val) { for(auto& kval: k) kval=val; }
  __host__ __device__ inline void operator+=(const TensorCoeffs& Tadd ) { for(int i=0;i<NCoeffs;i++) k[i]+= Tadd.k[i]; }
  __host__ __device__ inline void operator*=( const ftype mul ) {  for(auto& kval: k) kval*=mul; }
  inline __host__ __device__ ftype& get(const int ni, const int abc[]) {
    if(ni>Order) assert(0);
    const int a=abc[0], b=abc[1], c=abc[2], d=abc[3], e=abc[4];
    if(ni==0) return k[0]; else
    if(ni==1) return k[1 + a]; else
    if(ni==2) return k[1+DIM +a+b*DIM]; else
    if(ni==3) return k[1+DIM+DIM*DIM +a+b*DIM+c*DIM*DIM]; else
    if(ni==4) return k[1+DIM+DIM*DIM+DIM*DIM*DIM +a+b*DIM+c*DIM*DIM+d*DIM*DIM*DIM]; else
    if(ni==5) return k[1+DIM+DIM*DIM+DIM*DIM*DIM+DIM*DIM*DIM*DIM +a+b*DIM+c*DIM*DIM+d*DIM*DIM*DIM+e*DIM*DIM*DIM*DIM]; else
    assert(0);
    return k[0];
  }
};

inline __device__ void calc_moments_vec(const ftype4 gauge, const ftype f[LBMconsts::Qn], ftype mom_vec[LBMconsts::Qn] ) {
  for(int i=0; i<LBMconsts::Qn; i++) {
    mom_vec[i]=0;
    for(int j=0; j<LBMconsts::Qn; j++) {
      const int iq=j;
      const ftype3 v = LBMconsts::ef[iq]*gauge.w + make_ftype3(gauge.x,gauge.y,gauge.z);
      const ftype m_v = power(v.x,LBMconsts::MomentsPower[i].x) * power(v.y,LBMconsts::MomentsPower[i].y) * power(v.z,LBMconsts::MomentsPower[i].z);
      mom_vec[i]+= m_v*f[j];
    }
    //if(LBMconsts::MomentsPower[i].x+LBMconsts::MomentsPower[i].y+LBMconsts::MomentsPower[i].z>4) mom_vec[i]=0;
  }
}
template<int NORDER> inline __device__ void calc_moments_tensors(const ftype4 gauge, const ftype f[LBMconsts::Qn], TensorCoeffs<NORDER>& TC ) {
  for(int ni=0; ni<=TC.Order; ni++) {
    int Ncomps=1,Ndim=1, _ni=0; while( ++_ni <= ni ) { Ndim*=DIM; Ncomps+=Ndim; }
    for(int abc=0; abc<Ndim; abc++) {
      int abci[TC.Order>0?TC.Order:1]={0};
      for(int _i=0,Nds=1; _i<ni; _i++, Nds*=DIM ) abci[_i] = abc/Nds%DIM;
      ftype& an = TC.get(ni, abci);
      an=0;
      for(int iq=0; iq<LBMconsts::Qn; iq++) {
        const ftype3 v = LBMconsts::ef[iq]*gauge.w + make_ftype3(gauge.x,gauge.y,gauge.z);
        const ftype h_v = Hermite(ni, v, abci);
        an += h_v*f[iq];
      }
    }
  }
}
template<int NORDER> inline __device__ ftype eval_fi_Hermit(TensorCoeffs<NORDER>& TC, const int iq){
  ftype fi=0;
  for(int ni=0; ni<=TC.Order; ni++) {
    int Ncomps=1,Ndim=1, _ni=0, nfact=1; while( ++_ni <= ni ) { Ndim*=DIM; Ncomps+=Ndim; nfact*=_ni; }
    for(int abc=0; abc<Ndim; abc++) {
      int abci[TC.Order>0?TC.Order:1]={0};
      for(int _i=0,Nds=1; _i<ni; _i++, Nds*=DIM ) abci[_i] = abc/Nds%DIM;
      ftype& dn = TC.get(ni, abci);
      fi+= LBMconsts::w[iq]/nfact*dn*Hermite(ni, LBMconsts::ef[iq], abci);
    }
  }
  return fi;
}

template<int NORDER> inline __device__ TensorCoeffs<NORDER> convertAtoD(TensorCoeffs<NORDER>& an_tensor, const ftype4 gauge){
  const ftype sqT = gauge.w;
  const ftype dsqT = 1./sqT;
  const ftype T = sqT*sqT;
  const ftype u[3] = {gauge.x,gauge.y,gauge.z};
  TensorCoeffs<NORDER> dn_tensor;
  for(int ni=0; ni<=an_tensor.Order; ni++) {
    int Ncomps=1,Ndim=1, _ni=0; while( ++_ni <= ni ) { Ndim*=DIM; Ncomps+=Ndim; }
    for(int abc=0; abc<Ndim; abc++) {
      const int nO = an_tensor.Order>0?an_tensor.Order:1;
      int abci[MAX(nO,5)]={0};
      for(int _i=0,Nds=1; _i<ni; _i++, Nds*=DIM ) abci[_i] = abc/Nds%DIM;
      const int a=abci[0], b=abci[1], c=abci[2], d=abci[3], e=abci[4];
      ftype& dn  = dn_tensor.get(ni, abci);
      ftype& an  = an_tensor.get(ni, abci);
      int _0000[nO] = {0};
      /*int a000[nO]={a}  ;int b000[nO]={b}; int c000[nO]={c}; int d000[nO]={d}; int e000[nO]={e};
      int ab00[MAX(nO,2)]={a,b}; int ac00[MAX(nO,2)]={a,c}; int ad00[MAX(nO,2)]={a,d}; int bc00[MAX(nO,2)]={b,c}; int bd00[MAX(nO,2)]={b,d}; int cd00[MAX(nO,2)]={c,d};
      int bcd0[MAX(nO,3)]={b,c,d}; int acd0[MAX(nO,3)]={a,c,d}; int abd0[MAX(nO,3)]={a,b,d}; int abc0[MAX(nO,3)]={a,b,c};*/
      const ftype an0 = an_tensor.get(0, _0000);
      ftype An1_a,An1_b,An1_c,An1_d,An1_e,
            An2_ab,An2_ac,An2_ad,An2_ae,An2_bc,An2_bd,An2_be,An2_cd,An2_ce,An2_de,
            An3_abc,An3_abd,An3_abe,An3_acd,An3_ace,An3_ade,An3_bcd,An3_bce,An3_bde,An3_cde,
            An4_abcd,An4_abce,An4_abde,An4_acde,An4_bcde;
      if(ni>1) An1_a  = an_tensor.get(1,(int[]){a});
      if(ni>1) An1_b  = an_tensor.get(1,(int[]){b});
      if(ni>1) An1_c  = an_tensor.get(1,(int[]){c});
      if(ni>1) An1_d  = an_tensor.get(1,(int[]){d});
      if(ni>1) An1_e  = an_tensor.get(1,(int[]){e});
      if(NORDER>2)
      if(ni>2) { An2_ab = an_tensor.get(2,(int[]){a,b}); An2_ac = an_tensor.get(2,(int[]){a,c}); An2_ad = an_tensor.get(2,(int[]){a,d}); An2_ae = an_tensor.get(2,(int[]){a,e});
                 An2_bc = an_tensor.get(2,(int[]){b,c}); An2_bd = an_tensor.get(2,(int[]){b,d}); An2_be = an_tensor.get(2,(int[]){b,e});
                 An2_cd = an_tensor.get(2,(int[]){c,d}); An2_ce = an_tensor.get(2,(int[]){c,e}); An2_de = an_tensor.get(2,(int[]){d,e}); }
      if(NORDER>3)
      if(ni>3) { An3_abc = an_tensor.get(3,(int[]){a,b,c}); An3_abd = an_tensor.get(3,(int[]){a,b,d}); An3_abe = an_tensor.get(3,(int[]){a,b,e});
                  An3_acd = an_tensor.get(3,(int[]){a,c,d}); An3_ace = an_tensor.get(3,(int[]){a,c,e}); An3_ade = an_tensor.get(3,(int[]){a,d,e});
                  An3_bcd = an_tensor.get(3,(int[]){b,c,d}); An3_bce = an_tensor.get(3,(int[]){b,c,e}); An3_bde = an_tensor.get(3,(int[]){b,d,e});
                  An3_cde = an_tensor.get(3,(int[]){c,d,e}); }
      if(NORDER>4)
      if(ni>4) { An4_abcd = an_tensor.get(3,(int[]){a,b,c,d});
                 An4_abce = an_tensor.get(3,(int[]){a,b,c,e});
                 An4_abde = an_tensor.get(3,(int[]){a,b,d,e});
                 An4_acde = an_tensor.get(3,(int[]){a,c,d,e});
                 An4_bcde = an_tensor.get(3,(int[]){b,c,d,e}); }
      const ftype u1_ab = u[a]*u[b] - (T-1)*(a==b);
      const ftype u1_ac = u[a]*u[c] - (T-1)*(a==c);
      const ftype u1_ad = u[a]*u[d] - (T-1)*(a==d);
      const ftype u1_ae = u[a]*u[e] - (T-1)*(a==e);
      const ftype u1_bc = u[b]*u[c] - (T-1)*(b==c);
      const ftype u1_bd = u[b]*u[d] - (T-1)*(b==d);
      const ftype u1_be = u[b]*u[e] - (T-1)*(b==e);
      const ftype u1_cd = u[c]*u[d] - (T-1)*(c==d);
      const ftype u1_ce = u[c]*u[e] - (T-1)*(c==e);
      const ftype u1_de = u[d]*u[e] - (T-1)*(d==e);
      if(ni==0) dn = an; else
      if(ni==1) dn = dsqT*     ( an - u[a]*an0 ); else
      if(ni==2) dn = dsqT*dsqT*( an - u[a]*An1_b-u[b]*An1_a + (u[a]*u[b]-(T-1)*(a==b))*an0 ); else
      if(ni==3) dn = dsqT*dsqT*dsqT*( an - (u[a]*An2_bc+u[b]*An2_ac+u[c]*An2_ab) + (u1_ab*An1_c+u1_ac*An1_b+u1_bc*An1_a) 
                                         - ( u[a]*u[b]*u[c] - (T-1)*( u[a]*(b==c) + u[b]*(a==c) + u[c]*(a==b) ) )*an0 ); else
      if(NORDER>=4 && ni==4) dn = dsqT*dsqT*dsqT*dsqT*( an 
                           -( u[a]*An3_bcd + u[b]*An3_acd + u[c]*An3_abd + u[d]*An3_abc )
                           +( u1_ab*An2_cd + u1_ac*An2_bd + u1_ad*An2_bc + u1_bc*An2_ad + u1_bd*An2_ac + u1_cd*An2_ab )
                           -( u[a]*u[b]*u[c]*An1_d + u[a]*u[b]*u[d]*An1_c + u[a]*u[c]*u[d]*An1_b + u[b]*u[c]*u[d]*An1_a )
                           +(T-1)*( u[a]*(b==c)*An1_d + u[a]*(b==d)*An1_c + u[a]*(c==d)*An1_b
                                   +u[b]*(c==d)*An1_a + u[b]*(a==d)*An1_c + u[b]*(a==c)*An1_d
                                   +u[c]*(b==d)*An1_a + u[c]*(a==d)*An1_b + u[c]*(a==b)*An1_d
                                   +u[d]*(b==c)*An1_a + u[d]*(a==c)*An1_b + u[d]*(a==b)*An1_c )
                           + u[a]*u[b]*u[c]*u[d]*an0
                           -(T-1)*( u[a]*u[b]*(c==d) + u[a]*u[c]*(b==d) + u[a]*u[d]*(b==c) + u[b]*u[c]*(a==d) + u[b]*u[d]*(a==c) + u[c]*u[d]*(a==b) )*an0
                           +(T-1)*(T-1)*( (a==b)*(c==d) + (a==c)*(b==d) + (a==d)*(b==c) )*an0
                           ); else
      if(NORDER>=5 && ni==5) {
               dn = dsqT*dsqT*dsqT*dsqT*dsqT*( an 
                          -( u[a]*An4_bcde + u[b]*An4_acde + u[c]*An4_abde + u[d]*An4_abce + u[e]*An4_abcd )
                          +( u1_ab*An3_cde + u1_ac*An3_bde + u1_ad*An3_bce + u1_ae*An3_bcd + u1_bc*An3_ade + u1_bd*An3_ace + u1_be*An3_acd + u1_cd*An3_abe + u1_ce*An3_abd + u1_de*An3_abc )
                          -( u[a]*u[b]*u[c]*An2_de + u[a]*u[b]*u[d]*An2_ce + u[a]*u[b]*u[e]*An2_cd + u[a]*u[c]*u[d]*An2_be + u[a]*u[c]*u[e]*An2_bd + u[a]*u[d]*u[e]*An2_bc + u[b]*u[c]*u[d]*An2_ae + u[b]*u[c]*u[e]*An2_ad + u[b]*u[d]*u[e]*An2_ac + u[c]*u[d]*u[e]*An2_ab )
                    +(T-1)*(   u[a]*( (b==c)*An2_de + (b==d)*An2_ce + (b==e)*An2_cd + (c==d)*An2_be + (c==e)*An2_bd + (d==e)*An2_bc )
                             + u[b]*( (a==c)*An2_de + (a==d)*An2_ce + (a==e)*An2_cd + (c==d)*An2_ae + (c==e)*An2_ad + (d==e)*An2_ac )
                             + u[c]*( (a==b)*An2_de + (a==d)*An2_be + (a==e)*An2_bd + (b==d)*An2_ae + (b==e)*An2_ad + (d==e)*An2_ab )
                             + u[d]*( (a==b)*An2_ce + (a==c)*An2_be + (a==e)*An2_bc + (b==c)*An2_ae + (b==e)*An2_ac + (c==e)*An2_ab )
                             + u[e]*( (a==b)*An2_cd + (a==c)*An2_bd + (a==d)*An2_bc + (b==c)*An2_ad + (b==d)*An2_ac + (c==d)*An2_ab ) )
                          +( u[a]*u[b]*u[c]*u[d]*An1_e + u[a]*u[b]*u[c]*u[e]*An1_d + u[a]*u[b]*u[d]*u[e]*An1_c + u[a]*u[c]*u[d]*u[e]*An1_b + u[b]*u[c]*u[d]*u[e]*An1_a )
                    -(T-1)*( + u[a]*u[b]*( (c==d)*An1_e + (c==e)*An1_d + (d==e)*An1_c )
                             + u[a]*u[c]*( (b==d)*An1_e + (b==e)*An1_d + (d==e)*An1_b )
                             + u[a]*u[d]*( (b==c)*An1_e + (b==e)*An1_c + (c==e)*An1_b )
                             + u[a]*u[e]*( (b==c)*An1_d + (b==d)*An1_c + (c==d)*An1_b )
                             + u[b]*u[c]*( (a==d)*An1_e + (a==e)*An1_d + (e==d)*An1_a )
                             + u[b]*u[d]*( (a==c)*An1_e + (a==e)*An1_c + (e==c)*An1_a )
                             + u[b]*u[e]*( (a==c)*An1_d + (a==d)*An1_c + (d==c)*An1_a )
                             + u[c]*u[d]*( (a==b)*An1_e + (a==e)*An1_b + (b==e)*An1_a )
                             + u[c]*u[e]*( (a==b)*An1_d + (a==d)*An1_b + (b==d)*An1_a )
                             + u[d]*u[e]*( (a==b)*An1_c + (a==c)*An1_b + (b==c)*An1_a ) )
              +(T-1)*(T-1)*( + ( (a==b)*(c==d) + (a==c)*(b==d) + (a==d)*(b==c) )*An1_e
                             + ( (a==b)*(c==e) + (a==c)*(b==e) + (a==e)*(b==c) )*An1_d
                             + ( (a==b)*(d==e) + (a==d)*(b==e) + (a==e)*(b==d) )*An1_c
                             + ( (a==c)*(d==e) + (a==d)*(c==e) + (a==e)*(c==d) )*An1_b
                             + ( (b==c)*(d==e) + (b==d)*(c==e) + (b==e)*(c==d) )*An1_a )
                      ); } else
      //TODO fix this
      assert(0);
    }
  }
  return dn_tensor;
}

struct MomentsMatrix{
  ftype m[LBMconsts::Qn][LBMconsts::Qn*2];
  int InvStatus;

  __forceinline__ __device__ void init(const ftype4 gauge) {
    if(threadIdx.x!=0) return;
    for(int i=0; i<LBMconsts::Qn; i++) for(int j=0; j<LBMconsts::Qn; j++) {
      const int iq=j;
      ftype3 v = LBMconsts::ef[iq]*gauge.w + make_ftype3(gauge.x,gauge.y,gauge.z);
      m[i][j] = power(v.x,LBMconsts::MomentsPower[i].x) * power(v.y,LBMconsts::MomentsPower[i].y) * power(v.z,LBMconsts::MomentsPower[i].z);
    }
    InvStatus=0;
  }
  __forceinline__ __device__ ftype get(const int irow, const ftype f[LBMconsts::Qn]) {
    if(InvStatus==1) assert(0);
    ftype mom=0;
    for(int j=0; j<LBMconsts::Qn; j++) mom+= m[irow][j]*f[j];
    return mom;
  }
  __forceinline__ __device__ ftype get_inv(const int irow, const ftype m_vec[LBMconsts::Qn]) {
    if(InvStatus==0) assert(0);
    ftype fi=0;
    for(int j=0; j<LBMconsts::Qn; j++) {
      fi+= m[irow][j+LBMconsts::Qn]*m_vec[j];
    }
    return fi;
  }
  inline __device__ ftype operator* (const ftype f[LBMconsts::Qn]) { return 0; }

  inline __device__ void inverse() {
    if(threadIdx.x!=0) return;

    InvStatus^=1;

    #ifdef USE_DOUBLE
    const ftype err_tol=1e-13;
    #else
    const ftype err_tol=1e-30;
    #endif

    const int Qn = LBMconsts::Qn;
    //ftype m_t[Qn][2*Qn];
    for(int i=0;i<Qn;i++) {
      for(int j=0;j<2*Qn;j++) {
        //if(j<Qn) m_t[i][j] = m[i][j]; else
        if(j<Qn) {} else
        if(i==j-Qn) m[i][j]=1;
        else m[i][j]=0;
      }
    }
    for(int i=0;i<Qn;i++) {
      int temp = i;
      /* finding maximum i-th column element in last (dimension-j) rows */
      for (int j=i+1; j<Qn; j++) if(abs(m[j][i]) > abs(m[temp][i])) temp = j;
      if (abs(m[temp][i]) < err_tol) {
        printf("\n Elements are too small to deal with !!!( m[%d %d]=%g )\n ", temp,i,m[temp][i]);
        //        return -1;
      }
      /* swapping row which has maximum i-th column element */
      if(temp != i) for (int k=0; k<2*Qn; k++) {
        auto temporary = m[i][k];
        m[i][k] = m[temp][k];
        m[temp][k] = temporary;
      }
      for(int j=0;j<Qn;j++) {
        if(i!=j) {
          const ftype divm = m[i][i];
          ftype ratio = m[j][i]/divm;
          //if(divm<err_tol && divm>-err_tol) ratio=1/err_tol;
          for(int k=0; k<2*Qn; k++) m[j][k]-= ratio*m[i][k];
        } else {
          ftype divt = 1/m[i][i];
          for(int k=0; k<2*Qn; k++) m[i][k]*= divt;
        }
      }
    }
    //for(int i=0; i<Qn; i++) for(int j=0; j<Qn; j++) m[i][j] = m_t[i][Qn+j];
  }
};
