#ifndef IM2D_H
#define IM2D_H
#define GL_GLEXT_PROTOTYPES
#include <GL/glut.h>
#include <cuda_gl_interop.h>
#ifndef NOGL
#include <png.h>
#endif
#include "err.h"

extern int type_diag_flag;
//texture<float4, 1, cudaReadModeElementType> palTex; // 1D transfer function texture
void draw_scale();
struct image2D {
  GLuint bufferObj;
  cudaGraphicsResource* resource;
  int Nx, Ny;
  int xPos, yPos;
  uchar4* map4draw() {
    if(type_diag_flag>=3) printf("Отображаем буфер в адресное пространство\n");
    uchar4* devPtr; size_t size;
    if(CHECK_ERROR(cudaGraphicsMapResources(1, &resource, NULL))) throw(-1);
    if(CHECK_ERROR(cudaGraphicsResourceGetMappedPointer((void**) &devPtr, &size, resource))) throw(-1);
#ifdef COLOR_BACKGROUND
    uchar4* buf=new uchar4[Nx], col=make_uchar4(110,202,255,0);
    for(int i=0; i<Nx; i++) buf[i] = col;
    for(int i=0; i<Ny; i++)
  if(CHECK_ERROR(cudaMemcpy(devPtr+i*Nx, buf, 4*Nx, cudaMemcpyHostToDevice))) throw(-1);
    delete buf;
#else
    if(CHECK_ERROR(cudaMemset(devPtr, 200, size))) throw(-1);
#endif
    return devPtr;
  }
  void unmapAfterDraw() {
    if(type_diag_flag>=3) printf("Отключаем буфер от cuda для синхронизации с графикой\n");
    if(CHECK_ERROR(cudaGraphicsUnmapResources(1, &resource, NULL))) throw(-1);
    glutPostRedisplay();
  }
  void get_device(int pmajor=2, int pminor=1) {
    if(type_diag_flag>=1) printf("Находим подходящее устройство cuda и связываем его с GL устройством\n");
    cudaDeviceProp prop; memset(&prop, 0, sizeof(cudaDeviceProp));
    prop.major = pmajor; prop.minor = pminor;
    int dev;
    if(CHECK_ERROR(cudaChooseDevice(&dev, &prop))) throw(-1);
//    if(CHECK_ERROR(cudaGLSetGLDevice(dev))) throw(-1);
  }
  void init_image(int argc, char** argv, int _Nx, int _Ny, const char* tit) {
    Nx = _Nx; Ny = _Ny;
    if(type_diag_flag>=1); printf("Инициализация драйвера OpenGL\n");
    glutInit(&argc,argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowSize(Nx, Ny);
    //glutInitWindowPosition(100, 100);
    glutCreateWindow(tit);
    int wH=glutGet(GLUT_WINDOW_HEIGHT);
    yPos = wH<Ny?wH-Ny:0; xPos=0;

    if(type_diag_flag>=1) printf("Создаём пиксельный буфер в OpenGL, затем регистрируем его в cuda\n");
    glGenBuffers(1, &bufferObj);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, bufferObj);
    glBufferData(GL_PIXEL_UNPACK_BUFFER_ARB, Nx*Ny*sizeof(uchar4), NULL, GL_DYNAMIC_DRAW_ARB);
    if(CHECK_ERROR(cudaGraphicsGLRegisterBuffer(&resource, bufferObj, cudaGraphicsMapFlagsNone))) throw(-1);
  }
  void draw(char* tit) {
    glutSetWindowTitle(tit);
    glWindowPos2i(xPos,yPos);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, bufferObj);
    glDrawPixels(Nx,Ny, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
    draw_scale();
    glutSwapBuffers();
  }
  void clear() {
    if(type_diag_flag>=1) printf("Очистка OpenGL и CUDA\n");
    if(CHECK_ERROR(cudaGraphicsUnregisterResource(resource))) throw(-1);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
    glDeleteBuffers(1, &bufferObj);
  }
  ~image2D() {
  }
  void out2png(const char* pngName) {
#ifndef NOGL
    FILE* fp=fopen(pngName, "wb"); if(!fp) { printf("Не могу открыть файл %s\n", pngName); return; }
    uchar4* devPtr; size_t size;
    if(CHECK_ERROR(cudaGraphicsMapResources(1, &resource, NULL))) throw(-1);
    if(CHECK_ERROR(cudaGraphicsResourceGetMappedPointer((void**) &devPtr, &size, resource))) throw(-1);
    uchar4* hostPtr=(uchar4*)malloc(size);
    if(CHECK_ERROR(cudaMemcpy(hostPtr, devPtr, size, cudaMemcpyDeviceToHost))) throw(-1);
    
    glReadPixels(0, Ny-30, Nx, 30, GL_RGBA, GL_UNSIGNED_BYTE, hostPtr+Nx*(Ny-30));
    
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL,NULL);
    if(!png_ptr) { printf("png_create_write_struct failed\n"); return; }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr) { png_destroy_write_struct(&png_ptr, (png_infopp)NULL); return; }
    png_init_io(png_ptr, fp);
    png_set_filter(png_ptr, 0, PNG_FILTER_NONE);
    png_set_IHDR(png_ptr, info_ptr, Nx, Ny,
      8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_bytepp rows=(png_bytepp)malloc(sizeof(png_bytep)*Ny);
    uchar4* bmpp=hostPtr;
    for(int y=Ny-1; y>=0; y--){
      rows[y]=(png_bytep)malloc(sizeof(png_byte)*Nx*3);
      for(int x=0; x<Nx; x++, bmpp++) { rows[y][3*x] = bmpp->x; rows[y][3*x+1] = bmpp->y; rows[y][3*x+2] = bmpp->z; }
    }
    png_set_rows(png_ptr,info_ptr,rows);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    fclose(fp);
    
    for(int y=0; y<Ny; y++) free(rows[y]);
    free(rows);
    free(hostPtr);
    unmapAfterDraw();
    if(type_diag_flag>=0) printf("Картинка сохранена в %s\n", pngName);
#endif
  }
};
#endif//IM2D_H
