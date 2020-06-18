#ifndef FPAL_H
#define FPAL_H
#include <cuda.h>
#include <stdio.h>
#include "err.h"

// палитры творчески позаимствованы из aivlib
#define b 0.0001
#define q 0.25
#define h 0.75
#define H 0.5
const float fpal_array[] = { -1,
  1,0,H, h,0,0, 1,H,0, 1,1,0, 0,h,0, 0,0,0, 0,1,0, 0,0,1, H,0,1, 1,0,1, 0,1,1, -1,//hacked rainbow_pal
  //1,0,H, h,0,0, 1,H,0, 1,1,0, 0,h,0, 0,0,0, 0,0,1, 0,1,0, H,0,1, 1,0,1, 0,1,1, -1,//hacked rainbow_pal
  1,1,1, 1,0,H, h,0,0, 1,H,0, 1,1,0, 0,h,0, 0,1,1, 0,0,1, 1,0,1, H,0,1, 1,1,1, -1,//hacked rainbow_pal
  //H,H,H, 0,0,0, 1,0,0, 1,H,0, 1,1,0, 0,1,0, 0,1,1, 0,0,1, 1,0,1, 1,1,1, H,H,H, -1,//cyclic_pal
  1,0,H, h,0,0, 1,H,0, 1,1,0, 0,h,0, 0,0,1, H,0,1, 1,0,1, 0,1,1, 1,1,1, -1,//hacked rainbow_pal
  0,0,0, 1,0,0, 1,H,0, 1,1,0, 0,1,0, 0,1,1, 0,0,1, 1,0,1, 1,1,1, -1,//rainbow_pal
  1,1,1, 1,0,1, 0,0,1, 0,1,1, 0,1,0, 1,1,0, 1,H,0, 1,0,0, 0,0,0, -1,//inv_rainbow_pal
  q,0,0, 1,0,0, 1,H,0, 1,1,0, 1,1,1, 1,0,1, 0,0,1, 0,1,1, 0,1,0, -1,//neg_pos1_pal
  h,h,0, 1,1,0, 1,H,0, 1,0,0, 1,1,1, 0,1,0, 0,1,1, 0,0,1, 1,0,1, -1,//neg_pos2_pal
  1,1,1, 1,0,0, 1,H,0, 1,1,0, 0,1,0, 0,1,1, 0,0,1, 1,0,1, q,0,q, -1,//positive_pal
  1,1,1, 1,0,0, 1,H,0, 1,1,0, 0,1,0, 0,1,1, 0,0,1, 1,0,1, -1,//positive_pal
         1,0,0, 1,H,0, 1,1,0, 0,1,0, 0,1,1, 0,0,1, 1,0,1,        -1,//color_pal
         1,0,0, 1,1,0, 0,1,0, 0,1,1, 0,0,1, 1,0,1,        -1,//color_pal
         1,1,0, 1,0,0, 0,1,0, 0,0,1, 0,1,1,        -1,//color_pal
         1,0,0, 0,1,0, 0,0,1, 0,1,1,        -1,//color_pal
  0.4,1,1, 0.4,1,1, 0.4,1,1, 0.4,1,1, 0.28,0.7,1, 0,0,1,    b,b,b,   0,0.5,0, 0,1,0, 1,0,0, 1,0.5,0, 1,1,0, 1,1,1, -1, //zakirov_pal
  b,b,b, 1,1,1, 1,0,0, -1,//black_red_pal
  0,1,0, 1,1,1, 0,0,1, -1,//green_blue_pal
  1,1,1, 1,0,0, -1,//red_pal
  0.933,0.455,0, 1,1,1, 0,0,1, -1,//???pal
  //0,0,1, 0.9,0.9,0.9, 1,0,0, -1,//blue_red_pal
  0,0,1, 1,1,1, 1,0,0, -1,//blue_red_pal
  b,b,b, H,H,H, 1,1,1, -1,//grey_pal
  1,1,1, H,H,H, b,b,b, -1,//inv_grey_pal
  b,b,b, b,b,b, 1,1,1, -1,//black_black_white_pal
  1,0,H, h,0,0, 1,0,0, 1,H,0, 1,1,0, 0,h,0, 0,0,0, 0,1,0, 0,0,1, H,0,1, 1,0,1, 0,h,h, 0,1,1, -1,//super rainbow_pal
};
#undef H
#undef h
#undef q
#undef b
//прогрессия «в 10 раз за 10 кликов»:1.2   1.5     2     2.5     3      4      5       6      8     10
const float scale_step_array[] = { 6./5., 5./4., 4./3., 5./4., 6./5., 4./3., 5./4.,  6./5., 4./3., 5./4.};

//В нормальном коде так не делают, но в cuda текстуры должны иметь глобальные имена, так декларируем и соответствующие массивы так же
texture<float4, cudaTextureType1D, cudaReadModeElementType> fpal_col_tex;
texture<float, cudaTextureType1D> fpal_scale_tex;
cudaArray* fpal_col_texArray=0,* fpal_scale_texArray=0;

//------------------------------------
struct fpal_pars {
  float fmin, fmax, fscale, max_rgb;
  int start_pal, palID, palDim, pal3Daxis; float pscale;
  bool cyclic_pal, centric_pal, filter_pal, draw_flag, negate_flag, logFlag, transparency_discrete_flag;
  int scale_step, gamma_step, max_rgb_step, brightness_coff_step, transparency_mode, pal_scale_type;
  float gamma_pal, brightness_coff, base_val;
  //char ballast[128-76];//всего 128 байт
 public:
  void reset() {
    cyclic_pal = centric_pal = filter_pal = negate_flag = transparency_discrete_flag = false;
    start_pal = 1; palID = 0; palDim = 1; pal3Daxis = 2; pscale = 1.0f;
    transparency_mode = 1; pal_scale_type = 0;
    max_rgb = 1.0;
    draw_flag = false; logFlag = false;
    scale_step = gamma_step = max_rgb_step = 0; gamma_pal = 1.0;
    brightness_coff_step = 1; brightness_coff = scale_step_array[0]; base_val = 0.0;
  }
  void set_lim(float _fmin=0.0f, float _fmax=1.0f) {
    fmin = _fmin; fmax = _fmax;
    fscale = fmax>fmin?1.0f/(fmax-fmin):0.0;
  }
  void change_pal() {
    do { start_pal += 3; } while(fpal_array[start_pal] >= 0.0f);
    start_pal++; palID++;
    if(start_pal >= sizeof(fpal_array)/sizeof(float)) { start_pal = 1; palID = 0; }
  }
  void change_pal_back() {
    if(start_pal<=1) { start_pal = sizeof(fpal_array)/sizeof(float)-1; palID = 0; }
    else start_pal--;
    do { start_pal -= 3; } while(fpal_array[start_pal-1] >= 0.0f);
    palID--; //if(palID<0) palID = 0;
  }
  __device__ uchar4 invert_color(uchar4 c4) {
    return make_uchar4((c4.x+128)%256,(c4.y+128)%256,(c4.z+128)%256,255);
  }
  __device__ float4 get_color_norm_f4(float f, float br=1.0f) {
    f = 0.5f + pscale*f;
    return tex1D(fpal_col_tex,f);
  }
  __device__ float4 get_color_for3D(float f) {
  //const float pscale=im.pscale*0.01f, fscale=100.0f*im.fscale, fmin=0.5f-im.fmin*fscale;
    return tex1D(fpal_col_tex, 0.5f + pscale*0.01f*tex1D(fpal_scale_tex, 0.5f+(f-fmin)*100.0f*fscale));
  }
  __device__ float4 get_color_f4(float f) {
    float fn=(f-fmin)*fscale, br=1.0f;
    if(cyclic_pal && (fn<0.0f || fn>1.0f)) {
      float fi=floorf(fn);
      br = pow(brightness_coff,0.5f*fi);
      fn -= fi;
    }
    return get_color_norm_f4(0.01f*tex1D(fpal_scale_tex, 0.5f+100.0f*fn), br);
  }
  __device__ uchar4 get_color_norm(float f, float br=1.0) {
    f = 0.5f + pscale*f; br *= max_rgb;
    float4 col=tex1D(fpal_col_tex,f);
    return make_uchar4(__saturatef(br*col.x)*255, __saturatef(br*col.y)*255, __saturatef(br*col.z)*255, 255);
  }
  __device__ uchar4 get_color(float2 f) {
    float4 f4=get_color_for3D(f);
    return make_uchar4(__saturatef(f4.x)*255, __saturatef(f4.y)*255, __saturatef(f4.z)*255, 255);
  }
  __device__ float4 get_color_for3D(float2 f) {
    float phi=f.x, psi=f.y;
    float Cs=cos(phi), Sn=sin(phi);
    float cs=cos(psi), sn=sin(psi);
    return get_color_for3D(make_float4(Cs*cs,Sn*cs,sn, 0.0));
  }
  __device__ uchar4 get_color(float4 f) {
    float4 f4=get_color_for3D(f);
    return make_uchar4(__saturatef(f4.x)*255, __saturatef(f4.y)*255, __saturatef(f4.z)*255, __saturatef(f4.w)*255);
  }
  __device__ float4 get_color_for3D(float4 f) {
    for(int ia=0; ia<pal3Daxis; ia++) { float ft=f.x; f.x=f.y; f.y=f.z; f.z=ft; }
    //float3 v=make_float3(f.x*f.x,f.y*f.y,f.z*f.z);
    float A=length(make_float3(f.x,f.y,f.z));
    if(A==0) return make_float4(0.0f);
    float3 fs=make_float3(fabs(f.x),fabs(f.y),fabs(f.z)), v=fs; fs *= 1.0f/A;
    if(palID<-1) {
      if(fs.x<fs.y) { float t=fs.x; fs.x=fs.y; fs.y=t; }
      if(fs.y<fs.z) { float t=fs.y; fs.y=fs.z; fs.z=t; }
      if(fs.x<fs.y) { float t=fs.x; fs.x=fs.y; fs.y=t; }
      float cs=fs.x, sn=0.8+0.2*sqrt(fs.y*fs.y+fs.z*fs.z);
      for(int i=0; i<transparency_mode; i++) v *= cs;
      for(int i=0; i>transparency_mode; i--) v *= sn;
      float wv=1./(1+brightness_coff);
      const float3 r={0.0f,1.0f,1.0f}, g={1.0f,0.0f,1.0f}, b={1.0f,1.0f,0.0f}, w={wv,wv,wv};
      const float3 c={1.0f,0.0f,0.0f}, m={0.0f,1.0f,0.0f}, y={0.0f,0.0f,1.0f};
      float3 f3; int flag=-1-palID;
      if(flag&8) f3 = fscale*(v.x*(flag&1?w:(f.x<0?m:c))+v.y*(flag&2?w:(f.y<0?g:y))+v.z*(flag&4?w:(f.z<0?b:r)));
      else       f3 = fscale*(v.x*(flag&1?w:(f.x<0?c:m))+v.y*(flag&2?w:(f.y<0?y:g))+v.z*(flag&4?w:(f.z<0?r:b)));
      //float3 f3 = fscale*(v.x*(f.x<0?c:m)+v.y*(f.y<0?y:g)+v.z*(f.z<0?r:b));
      //float3 f3 = fscale*(v.x*(f.x<0?w:w)+v.y*(f.y<0?w:w)+v.z*(f.z<0?r:b));
      //float3 f3 = fscale*(fabs(f.x)*(f.x<0?R:r)+fabs(f.y)*(f.y<0?G:g)+fabs(f.z)*(f.z<0?B:b));
      return make_float4(f3.x,f3.y,f3.z, length(f3));
    } else if(palID>=0) {
      //return get_color_for3D(fabs(f.w));
      float3 f3, shadow=make_float3(tex1D(fpal_col_tex,0.5f*(pscale+1)+base_val));
      float cs=pow(v.z,4.0f/(1.0f+0.1f*abs(transparency_mode)));
      if(transparency_mode<=0 && f.z<0 || transparency_mode>=0 && f.z>0) f3 = cs*make_float3(get_color_for3D(f.w));
      else f3 = cs*shadow;
      f3 += (sqrt(v.x*v.x+v.y*v.y)/(1.0f+0.1f*brightness_coff))*shadow;
      return make_float4(f3.x,f3.y,f3.z, length(f3));
    } else {
      float phi=atan2(f.y,f.x), aphi=fabs(phi), z=f.z/A, L=z*pow(fabs(z),transparency_mode>=0?transparency_mode:1.0f/transparency_mode);
      //float psi=atan2(f.z,sqrt(f.y*f.y+f.x*f.x));
      float C=(1.0f-fabs(L))*fminf(1.0f,fscale*A), X=C*(1.0f-fabs(fmodf(phi+M_PI,2*M_PI/3)-M_PI/3)/(M_PI/3));
      float3 f3=make_float3(0);
      f3.x = aphi<M_PI/3.0?0.0f:(aphi>2.0*M_PI/3.0?C:X);
      f3.y = phi>M_PI/3.0?0.0f:(0>phi&&phi>-2.0*M_PI/3.0?C:X);
      f3.z = phi<-M_PI/3.0?0.0f:(0<phi&&phi<2.0*M_PI/3.0?C:X);
      float m=0.5*(L+1.0-C);
      return make_float4(m+f3.x,m+f3.y,m+f3.z, length(f3));
    }
  }
  __device__ uchar4 get_color(float f) {
    float fn=(f-fmin)*fscale, br=1.0f;
    if(cyclic_pal && (fn<0.0f || fn>1.0f)) {
      float fi=floorf(fn);
      br = pow(brightness_coff,0.5f*fi);
      fn -= fi;
    }
    return get_color_norm(0.01f*tex1D(fpal_scale_tex, 0.5f+100.0f*fn), br);
  }
  void bind2draw() {
    const int Nsc=100;
    float scale_data[Nsc+1], pal_data[128*4];
    const float fN=centric_pal?0.5f*Nsc:Nsc, dfN=1./fN;
    switch(pal_scale_type) {
     case 0: {//pow
      if(centric_pal) {
        for(int i=0; i<=Nsc/2; i++) {
          float v=fN*pow(i*dfN,gamma_pal);
          scale_data[Nsc/2+i] = fN+v;
          scale_data[Nsc/2-i] = fN-v;
        }
      } else { for(int i=0; i<=Nsc; i++) scale_data[i] = fN*pow(i*dfN,gamma_pal); }
     } break;
     case 1: {//exp
      float fNN=fN/(exp(gamma_pal)-1.0f);
      if(centric_pal) {
        for(int i=0; i<=Nsc/2; i++) {
          float v=fNN*(exp(i*dfN*gamma_pal)-1.0f);
          scale_data[Nsc-i] = 2*fN-v;
          scale_data[i] = v;
        }//for(int i=0; i<=Nsc; i++) printf("s[%d] = %g\n", i, double(scale_data[i]));
      } else { for(int i=0; i<=Nsc; i++) scale_data[i] = fNN*(exp(i*dfN*gamma_pal)-1.0f); }
     } break;
     default: for(int i=0; i<=Nsc; i++) scale_data[i] = i;
    };
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();
    if(CHECK_ERROR(cudaMallocArray(&fpal_scale_texArray, &channelDesc, Nsc+1))) throw(-1);
    if(CHECK_ERROR(cudaMemcpyToArray(fpal_scale_texArray, 0, 0, scale_data, (Nsc+1)*sizeof(float), cudaMemcpyHostToDevice))) throw(-1);
    fpal_scale_tex.normalized = false;
    fpal_scale_tex.filterMode =  cudaFilterModeLinear;
    fpal_scale_tex.addressMode[0] = cyclic_pal?cudaAddressModeBorder:cudaAddressModeClamp;
    if(CHECK_ERROR(cudaBindTextureToArray(fpal_scale_tex, fpal_scale_texArray))) throw(-1);
    int Nc; for(Nc=0; fpal_array[start_pal+3*Nc] >= 0.0f; Nc++);
    pscale = Nc-1.0f;
    for(int i=0; i<Nc; i++) {
      int ic=start_pal+3*i;
      if(!negate_flag || fpal_array[ic]+fpal_array[ic+1]+fpal_array[ic+2]==0.0) {
        pal_data[4*i] = fpal_array[ic]; pal_data[4*i+1] = fpal_array[ic+1]; pal_data[4*i+2] = fpal_array[ic+2]; pal_data[4*i+3] = 1.0f;
      } else {
        pal_data[4*i] = 1.0-fpal_array[ic]; pal_data[4*i+1] = 1.0-fpal_array[ic+1]; pal_data[4*i+2] = 1.0-fpal_array[ic+2]; pal_data[4*i+3] = 1.0f;
      }
    }
    //if(centric_pal) pal_data[4*(Nc/2)+3] = pal_data[4*((Nc-1)/2)+3] = 0.0f;
    //else pal_data[3] = pal_data[4*(Nc-1)+3] = 0.0f;
    if(transparency_discrete_flag) {
      for(int i=0; i<Nc; i++) pal_data[4*i+3] = (transparency_mode&(1<<i))==0?1.0f:0.0f;
    } else {
      float km=1.0*M_PI*transparency_mode, phi=0.0, inv=0.5;
      if(centric_pal) { km *= 2.0; phi = 0.5; }
      else if(transparency_mode<0) { phi = 1.0; inv = -inv; }
      if(transparency_mode>0) inv = -inv;
      for(int i=0; i<Nc; i++) pal_data[4*i+3] = pow(0.5+inv*cos(km*(i/(Nc-1.0)-phi)),brightness_coff);
    }
    cudaChannelFormatDesc channelDesc2 = cudaCreateChannelDesc<float4>();
    if(CHECK_ERROR(cudaMallocArray(&fpal_col_texArray, &channelDesc2, Nc))) throw(-1);
    if(CHECK_ERROR(cudaMemcpyToArray(fpal_col_texArray, 0, 0, pal_data, Nc*sizeof(float4), cudaMemcpyHostToDevice))) throw(-1);
    fpal_col_tex.normalized = false;//true;
    fpal_col_tex.filterMode = filter_pal?cudaFilterModePoint:cudaFilterModeLinear;
    fpal_col_tex.addressMode[0] = cudaAddressModeClamp;//cyclic_pal?cudaAddressModeWrap:cudaAddressModeClamp;
    if(CHECK_ERROR(cudaBindTextureToArray(fpal_col_tex, fpal_col_texArray, channelDesc2))) throw(-1);
  }
  void unbindAfterDraw() {
    if(CHECK_ERROR(cudaFreeArray(fpal_scale_texArray))) throw(-1);
    if(CHECK_ERROR(cudaFreeArray(fpal_col_texArray))) throw(-1);
  }
  void print_help() {
    printf("\
======= Палитра, цвета, масштабы:\n\
  p¦P \tВыбор палитры, в порядке уменьшения¦увеличения числа цветов (%.0f цвет[а/ов], начиная с поз. %d)\n\
«Ctl-p»\tПереключение вывода палитры в верх экрана (%d)\n\
  =¦- \tТочное (в 10 раз за 10 кликов) уменьшение¦увеличение пределов по цветовой оси (шаг %d)\n\
  +¦_ \tГрубое (в 10 раз) уменьшение¦увеличение пределов по цветовой оси (%g<f<%g)\n\
 0¦)¦(\tЦентрирование пределов палитры ¦ установка в 0 значения левого¦правого пределов\n\
  9¦8 \tУстановка пределов палитры в значения 1..9¦-1..1\n\
«Ctl-c»\tПереключение палитры из линейной в дискретную и обратно (%d). В частности - для уменьшения размера png\n\
  c¦C \tПереключение круговой¦центрированной палитры на ограниченную пределами и обратно (%d¦%d)\n\
  t¦T \tУвеличение¦уменьшение числа прозрачных цветов в палитре, (%d) д.б. делителем числа цветов в палитре (%d)\n\
«Ctl-t»\tПереключение на режим, в котором прозрачность каждого цвета в палитре кодируется битом числа %d (%d)\n\
  {¦} \tПереключение базового цвета в палитре (%g)\n\
«Ctl-g»\tПереключение функцию нелинейного шкалирования по цветовой оси степенная/показательная (%d) \n\
  [¦] \tУменьшение¦увеличение показателя степени при нелинейном шкалировании по цветовой оси (%g, шаг %d)\n\
  /¦\\\tУменьшение¦увеличение яркости цветов в палитре, важно для круговой (до %g, шаг %d)\n\
  ?¦| \tУменьшение¦увеличение относительной яркости разных циклов в круговой палитре (%g раз, шаг %d)\n\
«Ctr-r»\tСброс параметров в значения по умолчанию\n\
 «TAB»\tИнверсия цветов в палитре (кроме чёрного), для im3D решает проблему «грязи» на светлом фоне\n\
", pscale, start_pal, draw_flag, scale_step, fmin, fmax, filter_pal, cyclic_pal, centric_pal,
   transparency_mode, int(pscale)*(centric_pal?1:2), transparency_mode, transparency_discrete_flag,
   base_val, pal_scale_type, gamma_pal, gamma_step, max_rgb, max_rgb_step, pow(brightness_coff,0.5), brightness_coff_step);
  }
  bool key_func(unsigned char key, int x, int y) {
    //printf("key: %c\n", key);
    float fc=centric_pal?0.5*(fmax+fmin):0.0;
    switch(key) {
    case 18: reset(); break;
    case '-': set_lim(fc+(fmin-fc)*scale_step_array[scale_step], fc+(fmax-fc)*scale_step_array[scale_step]); scale_step = (scale_step+1)%10; break;
    case '=': scale_step = (scale_step+9)%10; set_lim(fc+(fmin-fc)/scale_step_array[scale_step], fc+(fmax-fc)/scale_step_array[scale_step]); break;
    case '_': set_lim(fc+(fmin-fc)*10.0f, fc+(fmax-fc)*10.0f); break;
    case '+': set_lim(fc+(fmin-fc)/10.0f, fc+(fmax-fc)/10.0f); break;
    case '8': centric_pal = true; set_lim(-1,1); break;
    case '9': set_lim(1,9); break;
    case '(': centric_pal = false; set_lim(fmin,0); break;
    case ')': centric_pal = false; set_lim(0,fmax); break;
    case '0': centric_pal = true; { float afmax=(fabs(fmin)>fabs(fmax))?fabs(fmin):fabs(fmax); set_lim(-afmax,afmax);} break;
    case 'c': cyclic_pal  ^= true; break;
    case 'C': centric_pal ^= true; break;
    case 't': transparency_mode++; break;
    case 'T': transparency_mode--; break;
    case  20: transparency_mode = 0; transparency_discrete_flag ^= true; break;
    case 3: filter_pal ^= true; break;
    case 'p': change_pal(); break;
    case 'P': change_pal_back(); break;
    case  16: draw_flag ^= true; break;
    case  9: negate_flag ^= true; break;
    case  7: pal_scale_type = (pal_scale_type+1)%2; break;
    case '\\':max_rgb *= sqrt(scale_step_array[max_rgb_step]); max_rgb_step = (max_rgb_step+1)%10; break;
    case '/': max_rgb_step = (max_rgb_step+9)%10; max_rgb /= sqrt(scale_step_array[max_rgb_step]); break;
    case '|': brightness_coff *= scale_step_array[brightness_coff_step]; brightness_coff_step = (brightness_coff_step+1)%10; break;
    case '?': brightness_coff_step = (brightness_coff_step+9)%10; brightness_coff /= scale_step_array[brightness_coff_step]; break;
    case 'l': logFlag=~logFlag; break;
    case 'L': logFlag=false; break;
    case '[': gamma_pal *= scale_step_array[gamma_step]; gamma_step = (gamma_step+1)%10; break;
    case ']': gamma_step = (gamma_step+9)%10; gamma_pal /= scale_step_array[gamma_step]; break;
    case '{': base_val -= 0.5; break;
    case '}': base_val += 0.5; break;
    default : return false;
    } return true;
  }
};

//------------------------------------
struct image_pars: public fpal_pars {
  uchar4* bmp,* bmp4backgrownd;
  unsigned int nFrame;
  unsigned int nFpng;
  void reset() {
    fpal_pars::reset();
    nFrame = 0;
    nFpng = 0;
    bmp = bmp4backgrownd = 0;
  }
};
#endif//FPAL_H
