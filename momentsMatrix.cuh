inline __device__ void calc_moments_vec(const ftype4 gauge, const ftype f[LBMconsts::Qn], ftype mom_vec[LBMconsts::Qn] ) {
  for(int i=0; i<LBMconsts::Qn; i++) {
    mom_vec[i]=0;
    for(int j=0; j<LBMconsts::Qn; j++) {
      const int iq=j;
      const ftype3 v = LBMconsts::ef[iq]*gauge.w + make_ftype3(gauge.x,gauge.y,gauge.z);
      const ftype m_v = power(v.x,LBMconsts::MomentsPower[i].x) * power(v.y,LBMconsts::MomentsPower[i].y) * power(v.z,LBMconsts::MomentsPower[i].z);
      mom_vec[i]+= m_v*f[j];
    }
  }
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
