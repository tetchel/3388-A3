
/*            PURPOSE : Simple framework for ray-tracing

        PREREQUISITES : matrix.h
 */
#define _GNU_SOURCE
#include <X11/Xlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#define SPHERE   1
#define PLANE    2

#define EPSILON 0.00001
#define MAX_OBJECTS 100
#define MAX_INTENSITY 255.0

#define Ex 10.0 /* position of camera frame of reference */
#define Ey 10.0
#define Ez 10.0

#define Gx 0.0 /* camera gaze direction */
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define Lx 10.0 /* light position */
#define Ly -10.0
#define Lz 10.0

#define Near 1.0
#define Far 25.0

#define THETA 45.0
#define ASPECT 1.5

#define H 800 /* window height in pixels */

typedef struct {
    int width, height ;
} window_t ;

typedef struct {
    dmatrix_t UP ;
    dmatrix_t E ;
    dmatrix_t G ;
    dmatrix_t u, v, n ;
} camera_t ;

typedef struct {
    double r, g, b ;
} color_t ;

typedef struct {
    int type ;
    double (*intersection_function)(dmatrix_t *,dmatrix_t *) ;
    dmatrix_t M, Minv ;
    color_t specular_color, diffuse_color, ambient_color ;
    double reflectivity, specular_coeff, diffuse_coeff, ambient_coeff, f ;
} object_t ;

typedef struct {
    dmatrix_t position ;
    color_t color ;
    color_t intensity ;
} light_t ;

object_t object[MAX_OBJECTS] ;
int nobjects = 0 ;


Display *InitX(Display *d, Window *w, int *s, window_t *Window) {

    d = XOpenDisplay(NULL) ;
    if(d == NULL) {
        printf("Cannot open display\n") ;
        exit(1) ;
    }
    *s = DefaultScreen(d) ;
    *w = XCreateSimpleWindow(d,RootWindow(d,*s),0,0,Window->width,Window->height,1,BlackPixel(d,*s),WhitePixel(d, *s)) ;
    Atom delWindow = XInternAtom(d,"WM_DELETE_WINDOW",0) ;
    XSetWMProtocols(d,*w,&delWindow,1) ;
    XSelectInput(d,*w,ExposureMask | KeyPressMask) ;
    XMapWindow(d,*w) ;
    return(d) ;
}

void SetCurrentColorX(Display *d, GC *gc, unsigned int r, unsigned int g, unsigned int b) {

    XSetForeground(d,*gc,r << 16 | g << 8 | b) ;
}

void SetPixelX(Display *d, Window w, int s, int i, int j) {

    XDrawPoint(d,w,DefaultGC(d,s),i,j) ;
}


void QuitX(Display *d, Window w) {

    XDestroyWindow(d,w) ;
    XCloseDisplay(d) ;
}

double dot_product(dmatrix_t *a, dmatrix_t *b) {

    double s  ;

    s = a->m[1][1]*b->m[1][1] + a->m[2][1]*b->m[2][1] + a->m[3][1]*b->m[3][1] ;

    return(s) ;
}

light_t *build_light(light_t *light, dmatrix_t *position, color_t color, color_t intensity) {

    dmat_alloc(&light->position,4,1) ;
    light->position = *position ;
    light->color.r = color.r ;
    light->color.g = color.g ;
    light->color.b = color.b ;
    light->intensity.r = intensity.r ;
    light->intensity.g = intensity.g ;
    light->intensity.b = intensity.b ;
    return light ;
}

window_t *build_window(window_t *Window, int height, float aspect) {

    Window->height = height ;
    Window->width =  aspect*height ;

    return(Window) ;
}

camera_t *build_camera(camera_t *Camera, window_t *Window) {

    dmat_alloc(&Camera->E,4,1) ;

    Camera->E.m[1][1] = Ex ;
    Camera->E.m[2][1] = Ey ;
    Camera->E.m[3][1] = Ez ;
    Camera->E.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->G,4,1) ;

    Camera->G.m[1][1] = Gx ;
    Camera->G.m[2][1] = Gy ;
    Camera->G.m[3][1] = Gz ;
    Camera->G.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->n,4,1) ;
    Camera->n = *dmat_normalize(dmat_sub(&Camera->E,&Camera->G)) ;
    Camera->n.l = 3 ;

    dmat_alloc(&Camera->UP,4,1) ;

    Camera->UP.l = 3 ;

    Camera->UP.m[1][1] = UPx ;
    Camera->UP.m[2][1] = UPy ;
    Camera->UP.m[3][1] = UPz ;
    Camera->UP.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->u,4,1) ;

    Camera->u = *dmat_normalize(dcross_product(&Camera->UP,&Camera->n)) ;
    Camera->v = *dmat_normalize(dcross_product(&Camera->n,&Camera->u)) ;

    return(Camera) ;
}

double sphere_intersection(dmatrix_t *e, dmatrix_t *d ) {

    return 0.0 ;
}

dmatrix_t *intersection_coordinates(dmatrix_t *e, dmatrix_t *direction, double t) {

    /* returns the 3D intersection point, given a t value for an intersection of an object with the ray */

    dmatrix_t *intersection ;

    intersection = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(intersection,4,1) ;

    /* intersection->m[1][1] = e->m[1][1] + direction->m[1][1]*t ;
    intersection->m[2][1] = e->m[2][1] + direction->m[2][1]*t ;
    intersection->m[3][1] = e->m[3][1] + direction->m[3][1]*t ;
    intersection->m[4][1] = 1.0 ; */

    return intersection ;
}

dmatrix_t *normal_to_surface(object_t *object, dmatrix_t *intersection) {

    /* returns normal surface unit vector at intersection point */

    dmatrix_t *normal ;

    normal = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(normal,4,1) ;

    /* if ((*object).type == 1) {
        normal->m[1][1] = intersection->m[1][1] ;
        normal->m[2][1] = intersection->m[2][1] ;
        normal->m[3][1] = intersection->m[3][1] ;
        normal->m[4][1] = 0.0 ;
    } */
    return normal ;
}

int find_min_hit_time(double t0[N_OBJECTS]) {

    /* finds the smallest t value for an intersection with an object */

    return 0 ;
}

dmatrix_t *ray_direction(camera_t *Camera, window_t *Window, double height, double width, int i, int j) {

    dmatrix_t *d ;

    d = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(d,4,1) ;

    d->m[1][1] = -1.0*Near*Camera->n.m[1][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[1][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[1][1] ;

    d->m[2][1] = -1.0*Near*Camera->n.m[2][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[2][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[2][1] ;

    d->m[3][1] = -1.0*Near*Camera->n.m[3][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[3][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[3][1] ;

    d->m[4][1] = 0.0 ;

    return(d) ;
}

dmatrix_t *vector_to_light_source(dmatrix_t *intersection, dmatrix_t *light_position) {

    /* returns a unit vector to light source at intersection point */

    dmatrix_t *s ;

    s = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(s,4,1) ;

    return s ;
}

dmatrix_t *vector_to_center_of_projection(dmatrix_t *intersection, dmatrix_t *e) {

    /* returns a unit vector to origin of camera frame of reference at intersection point */

    dmatrix_t *v ;

    v = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(v,4,1) ;

    return v ;
}

dmatrix_t *vector_to_specular_reflection(dmatrix_t *N, dmatrix_t *S) {

    /* returns a unit vector in the direction of the specular reflection given a surface normal and a light source */

    dmatrix_t *r ;

    r = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(r,4,1) ;

    return r ;

}

color_t color_init(double r, double g, double b) {

    color_t s ;

    s.r = r ;
    s.g = g ;
    s.b = b ;

    return s ;
}

color_t color_mult(double a, color_t c) {

    color_t s ;

    s.r = a*c.r ;
    s.g = a*c.g ;
    s.b = a*c.b ;

    return s ;
}

color_t color_add(color_t c1, color_t c2) {

    color_t s ;

    s.r = c1.r + c2.r ;
    s.g = c1.g + c2.g ;
    s.b = c1.b + c2.b ;

    return s ;
}

color_t shade(light_t *light,       //light object
              object_t *object,     //array of ALL objects
              dmatrix_t *e,         //the origin of viewing system expressed in world system
              dmatrix_t *d,         //vector representing the current ray
              color_t color,        //no idea what this is for ("old" colour?)
              color_t background,   //colour of the background (why is this needed?)
              int level,            //"closeness" to the near plane for occlusion - does higher level mean closer? seems to
              int i, int j)         //coordinates of pixel (i,j)
              {

    //invocation: pixel = shade(&light,object,&Camera.E,&direction,pixel,background,2,i,j) ;

    /* main ray-tracing routine. given a ray, performs the following:

        for all objects in the scene {
            transforms the ray and the camera eye position with Minv of the object
            collects t values of ray intersection with generic object
        }*/
        //what is minv
        for(int i=0; i < nobjects; i++) {

        }

    /*
        if one or more intersections found {
            finds the intersection with minimum t value
            transforms the light with Minv
            computes the coordinates of the intersection
            using the coordinates of the intersection {
                computes the unit surface normal vector
                computes unit vector to center of the camera
                computes unit vector to light source
                computes vector of specular reflection
            }
            computes specular, diffuse, and ambient light components
            returns the R,G,B intensity of the pixel
     } */

    color.r = 0.0 ;
    color.g = 0.0 ;
    color.b = 0.0 ;

    return color ;
}

object_t *build_object(int object_type, dmatrix_t *M, color_t ambient_color, color_t diffuse_color, color_t specular_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double f, double reflectivity) {

    object_t object ;

    object.M = *M ;
    dmat_alloc(&object.Minv,4,4) ;
    object.Minv = *dmat_inverse(&object.M) ;

    object.reflectivity = reflectivity ;

    object.specular_color.r = specular_color.r ;
    object.specular_color.g = specular_color.g ;
    object.specular_color.b = specular_color.b ;
    object.specular_coeff = specular_coeff ;
    object.f = f ;

    object.diffuse_color.r = diffuse_color.r ;
    object.diffuse_color.g = diffuse_color.g ;
    object.diffuse_color.b = diffuse_color.b ;
    object.diffuse_coeff = diffuse_coeff ;

    object.ambient_color.r = ambient_color.r ;
    object.ambient_color.g = ambient_color.g ;
    object.ambient_color.b = ambient_color.b ;
    object.ambient_coeff = ambient_coeff ;

    switch (object_type) {

        case SPHERE :   object.type = SPHERE ;
                        object.intersection_function = &sphere_intersection ;
            break ;
    }
    nobjects++ ;
    //is this ok? will it not be destroyed after return?
    //object should be malloced
    return(&object) ;

}

int main() {

    Display *d ;
    Window w ;
    XEvent e ;

    int i, j, s ;

    camera_t Camera ;
    window_t Window ;
    light_t light ;
    dmatrix_t M, light_position ;
    color_t pixel, background, light_intensity, light_color, ambient_color, diffuse_color, specular_color ;
    double height, width, aspect, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity ;

    /* Set the background color */

    background.r = 0.0 ;
    background.g = 0.0 ;
    background.b = 0.0 ;

    /* Set up light position, intensity, and color */

    dmat_alloc(&light_position,4,1) ;

    light_position.m[1][1] = Lx ;
    light_position.m[2][1] = Ly ;
    light_position.m[3][1] = Lz ;
    light_position.m[4][1] = 1.0 ;

    light_intensity.r = 1.0 ;
    light_intensity.g = 1.0 ;
    light_intensity.b = 1.0 ;

    light_color.r = 1.0 ;
    light_color.g = 1.0 ;
    light_color.b = 1.0 ;

    light = *build_light(&light,&light_position,light_color,light_intensity) ;

    /* Build display window and synthetic camera */

    Window = *build_window(&Window,H,ASPECT) ;
    Camera = *build_camera(&Camera,&Window) ;

    /* Create a first sphere */

     dmat_alloc(&M,4,4) ;
    M = *dmat_identity(&M) ;

    M.m[1][4] = 0.0 ; /* sphere at center of the world coordinates */
    M.m[2][4] = 0.0 ;

    reflectivity = 0.2 ;

    specular_color.r = 1.0 ;
    specular_color.g = 1.0 ;
    specular_color.b = 1.0 ;
    specular_coeff = 0.4 ;
    f = 10.0 ;

    diffuse_color.r = 0.0 ;
    diffuse_color.g = 0.0 ;
    diffuse_color.b = 1.0 ;
    diffuse_coeff = 0.4 ;

    ambient_color.r = 0.0 ;
    ambient_color.g = 0.0 ;
    ambient_color.b = 1.0 ;
    ambient_coeff = 0.2 ;

    object[nobjects] = *build_object(SPHERE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;

    /* Set near plane dimensions */

    aspect = ASPECT ;
    height = Near*tan(M_PI/180.0 * THETA/2.0) ;
    width = height*aspect ;

    dmatrix_t direction ;
    dmat_alloc(&direction,4,1) ;

    d = InitX(d,&w,&s,&Window) ;
    XNextEvent(d, &e) ;

    while (1) {
        XNextEvent(d, &e) ;
        if (e.type == Expose) {

            /* ray tracer main loop */

            for (i = 0 ; i < Window.width ; i++) {
                for (j = 0  ; j < Window.height ; j++) {
                    direction = *ray_direction(&Camera,&Window,height,width,i,j) ;
                    pixel = shade(&light,object,&Camera.E,&direction,pixel,background,2,i,j) ;
                    SetCurrentColorX(d,&(DefaultGC(d,s)),(int)pixel.r,(int)pixel.g,(int)pixel.b) ;
                    SetPixelX(d,w,s,i,Window.height - (j + 1)) ;
                }
            }
        }
        if (e.type == KeyPress)
            break ;
        if (e.type == ClientMessage)
            break ;
    }
    free_dmatrix(direction.m,1,direction.l,1,direction.c) ;
    QuitX(d,w) ;
}
