/*
 CSCI 480
 Assignment 3 Raytracer
 
 Name: <Your name here>
 */

#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include "pic.h"
#include <string.h>
#include <math.h>
#include <iostream>
using namespace std;

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct point {
    double x;
    double y;
    double z;
};

/* Class to represent A single vector's values and Standard Vector Functions*/
struct vector: public point{
    
    static vector cross_product (vector v1, vector v2) {
        vector v3;
        v3.x = v1.y*v2.z - v1.z*v2.y;
        v3.y = v2.x*v1.z - v1.x*v2.z;
        v3.z = v1.x*v2.y - v1.y*v2.x;
        return v3;
    }
    
    static float dot_product (point v1, point v2){
        return v1.x * v2.x +  v1.y * v2.y +  v1.z * v2.z;
    }
    
    void normalize (){
        double vectorLength = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        this->x = this->x/vectorLength;
        this->y = this->y/vectorLength;
        this->z = this->z/vectorLength;
    }
    
    static vector determine_vector (point p1, point p2){
        vector v;
        v.x = p2.x-p1.x;
        v.y = p2.y-p1.y;
        v.z = p2.z-p1.z;
        v.normalize();
        return v;
    }
};

struct Vertex
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double normal[3];
    double shininess;
};

typedef struct _Triangle
{
    struct Vertex v[3];
    
    float compute_area (point p1, point p2, point p3){
        float a = sqrt(p1.x * p2.x + p1.y * p2.y + p1.z * p2.z);
        float b = sqrt(p3.x * p2.x + p3.y * p2.y + p3.z * p2.z);
        float c = sqrt(p1.x * p3.x + p1.y * p3.y + p1.z * p3.z);
        float s = (a + b + c)/2;
        return sqrt(s * (s-a) * (s-b) * (s-c));
    }
} Triangle;

typedef struct _Sphere
{
    double position[3];
    double color_diffuse[3];
    double color_specular[3];
    double shininess;
    double radius;
} Sphere;

typedef struct _Light
{
    double position[3];
    double color[3];
} Light;

struct Ray
{
    point origin;
    vector direction;
    
    Ray (point p, vector v){
        origin = p; direction = v;
    }
    
    float check_sphere_intersection (Sphere sphere1){
        //(x0 + xd t - xc)2 + (y0 + yd t - yc)2 + (z0 + zd t - zc)2 - r2 = 0
        //Simplifies to at2 + bt + c = 0
            //a = xd2 + yd2+ zd2 = 1
            //b = 2(xd(x0-xc)+yd(y0-yc)+zd(z0-zc))
            //c = (x0-xc)2 + (y0-yc)2 +(z0-zc)2 - r2
        
        float b, c;
        b = 2*((direction.x * (origin.x - sphere1.position[0])) +
            (direction.y * (origin.y - sphere1.position[1])) +
            (direction.z * (origin.z - sphere1.position[2])));
        c = pow(origin.x-sphere1.position[0], 2) + pow(origin.y-sphere1.position[1], 2) + pow(origin.z-sphere1.position[2], 2) - pow(sphere1.radius, 2);
        
        //Solve to get t0, t1
        //t0,1 = (-b +/- sqrt(b2-4c))/2
        //If t0, t1 > 0, return min (t0, t1)
        
        float t0, t1;
        t0 = (-b + sqrt(pow(b, 2) - 4*c))/2;
        t1 = (-b - sqrt(pow(b, 2) - 4*c))/2;
        if (t0 < 0 && t1 < 0)
            return -1;
        else
            return fmax(t0, t1);
    }
    
    float check_triangle_intersection (Triangle triangle1){
        //Check if ray intersects triangle plane
        
        point p1, p2, p3;
        p1.x = triangle1.v[0].position[0];
        p1.y = triangle1.v[0].position[1];
        p1.z = triangle1.v[0].position[2];
        p2.x = triangle1.v[1].position[0];
        p2.y = triangle1.v[1].position[1];
        p2.z = triangle1.v[1].position[2];
        p3.x = triangle1.v[2].position[0];
        p3.y = triangle1.v[2].position[1];
        p3.z = triangle1.v[2].position[2];
        
        vector A, B;
        A.x = p2.x - p1.x;
        A.y = p2.y - p1.y;
        A.z = p2.z - p1.z;
        B.x = p3.x - p1.x;
        B.y = p3.y - p1.y;
        B.z = p3.z - p1.z;
        
        //Normal = A x B
        vector normal = vector::cross_product(A, B);
        vector d;
        
        /*
         
         t = n . (v0 - p0) / n . (p1-p0)
         
         point t = -(ax0 + by0 + cz0 + d)/ axd + byd + czd =
         -(n . p0 + d) / n . d
         
         if n . d = 0, no intersection
         if t <= 0 intersection is behind ray origin
         */
        
        point intersection;
        intersection.x = p1.x - origin.x;
        intersection.y = p1.y - origin.y;
        intersection.z = p1.z - origin.z;
        
        float t = vector::dot_product(normal, intersection) / (vector::dot_product(normal, direction));
        
        intersection.x = p1.x + t*direction.x;
        intersection.y = p1.x + t*direction.y;
        intersection.z = p1.x + t*direction.z;
        
        /*
         Project the point and triangle onto a plane
         • Pick a plane not perpendicular to triangle (such
         a choice always exists)
         • x = 0, y = 0, or z = 0
         2. Then, do the 2D test in the plane, by computing
         barycentric coordinates (follows next)
         
         Now we have 3 points instead of 2
         • Define 3 barycentric coordinates α, β, γ
         • p = α p1 + β p2 + γ p3
         • p inside triangle iff 0 ≤ α, β, γ ≤ 1, α + β + γ = 1
         • How do we calculate α, β, γ?
         
         Coordinates are ratios of triangle areas
         • Areas in these formulas should be signed
         - Clockwise (-) or anti-clockwise (+) orientation of the triangle
         - Important for point-in-triangle test
         */
        
        float alpha,beta,gamma;
       

        alpha = triangle1.compute_area(intersection, p2, p3) / triangle1.compute_area(p1, p2, p3);
        beta = triangle1.compute_area(p1, intersection, p3) / triangle1.compute_area(p1, p2, p3);
        gamma = triangle1.compute_area(p1, p2, intersection) / triangle1.compute_area(p1, p2, p3);
        if (0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && 0 <= gamma && gamma <= 1)
            return 1;
        else
            return -1;
    }
    
};
Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//MODIFY THIS FUNCTION
void draw_scene()
{
    unsigned int x,y;
    //simple output
    for(x=0; x<WIDTH; x++)
    {
        glPointSize(2.0);
        glBegin(GL_POINTS);
        for(y=0;y < HEIGHT;y++)
        {
            plot_pixel(x,y,x%256,y%256,(x+y)%256);
        }
        glEnd();
        glFlush();
    }
    printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
    glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
    glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
    buffer[HEIGHT-y-1][x][0]=r;
    buffer[HEIGHT-y-1][x][1]=g;
    buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
    point camera;
    camera.x = 0; camera.y = 0; camera.z = 0;
    point screen;
    screen.x = x; screen.y = y; screen.z = -2;
    point screen2 = screen;
    screen2.z = -3;
    vector direction = vector::determine_vector(screen, screen2);
    Ray r1(screen, direction);
    

    bool intersects = false;
    //For every object in the scene
    //Check spheres
    for (int i = 0; i < num_spheres; i++){
        //If they intersect
        if (r1.check_sphere_intersection(spheres[i]) >= 0){
            //Plot Object's Pixel
            intersects = true;
        }
    }
    //Check Triangles
    for (int i = 0; i < num_triangles; i++){
        //If they intersect
        if (r1.check_triangle_intersection(triangles[i]) >= 0){
            //Plot Object's Pixel
            intersects = true;
        }
    }
    if (intersects)
         plot_pixel_display(x,y,256.0,0.0,0.0);
    else
        plot_pixel_display(x,y,0.0,0.0,0.0);
    
    
   // if(mode == MODE_JPEG)
     //   plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
    Pic *in = NULL;
    
    in = pic_alloc(640, 480, 3, NULL);
    printf("Saving JPEG file: %s\n", filename);
    
    memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
    if (jpeg_write(filename, in))
        printf("File saved Successfully\n");
    else
        printf("Error in Saving\n");
    
    pic_free(in);
    
}

void parse_check(char *expected,char *found)
{
    if(strcasecmp(expected,found))
    {
        char error[100];
        printf("Expected '%s ' found '%s '\n",expected,found);
        printf("Parse error, abnormal abortion\n");
        exit(0);
    }
    
}

void parse_doubles(FILE*file, char *check, double p[3])
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check(check,str);
    fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
    printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
    char str[100];
    fscanf(file,"%s",str);
    parse_check("rad:",str);
    fscanf(file,"%lf",r);
    printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
    char s[100];
    fscanf(file,"%s",s);
    parse_check("shi:",s);
    fscanf(file,"%lf",shi);
    printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
    argv = "/Users/joshgreenberger/Documents/Class Material/Senior/Graphics/assign3/RayTracer/screenfile.txt";
    FILE *file = fopen(argv,"r");
    if (file == NULL) {
        printf ("Can't open file.\n");
        exit(1);
    }
    int number_of_objects;
    char type[50];
    int i;
    Triangle t;
    Sphere s;
    Light l;
    fscanf(file,"%i",&number_of_objects);
    
    printf("number of objects: %i\n",number_of_objects);
    char str[200];
    
    parse_doubles(file,"amb:",ambient_light);
    
    for(i=0;i < number_of_objects;i++)
    {
        fscanf(file,"%s\n",type);
        printf("%s\n",type);
        if(strcasecmp(type,"triangle")==0)
        {
            
            printf("found triangle\n");
            int j;
            
            for(j=0;j < 3;j++)
            {
                parse_doubles(file,"pos:",t.v[j].position);
                parse_doubles(file,"nor:",t.v[j].normal);
                parse_doubles(file,"dif:",t.v[j].color_diffuse);
                parse_doubles(file,"spe:",t.v[j].color_specular);
                parse_shi(file,&t.v[j].shininess);
            }
            
            if(num_triangles == MAX_TRIANGLES)
            {
                printf("too many triangles, you should increase MAX_TRIANGLES!\n");
                exit(0);
            }
            triangles[num_triangles++] = t;
        }
        else if(strcasecmp(type,"sphere")==0)
        {
            printf("found sphere\n");
            
            parse_doubles(file,"pos:",s.position);
            parse_rad(file,&s.radius);
            parse_doubles(file,"dif:",s.color_diffuse);
            parse_doubles(file,"spe:",s.color_specular);
            parse_shi(file,&s.shininess);
            
            if(num_spheres == MAX_SPHERES)
            {
                printf("too many spheres, you should increase MAX_SPHERES!\n");
                exit(0);
            }
            spheres[num_spheres++] = s;
            spheres[num_spheres-1].position[0] += WIDTH/2;  //set to center of screen
            spheres[num_spheres-1].position[1] += HEIGHT/2;
            spheres[num_spheres-1].radius *= 50;        //set size to screen
        }
        else if(strcasecmp(type,"light")==0)
        {
            printf("found light\n");
            parse_doubles(file,"pos:",l.position);
            parse_doubles(file,"col:",l.color);
            
            if(num_lights == MAX_LIGHTS)
            {
                printf("too many lights, you should increase MAX_LIGHTS!\n");
                exit(0);
            }
            lights[num_lights++] = l;
        }
        else
        {
            printf("unknown type in scene description:\n%s\n",type);
            exit(0);
        }
    }
    return 0;
}

void display()
{
    
}

void init()
{
    glMatrixMode(GL_PROJECTION);
    glOrtho(0,WIDTH,0,HEIGHT,1,-1);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
    //hack to make it only draw once
    static int once=0;
    if(!once)
    {
        draw_scene();
        if(mode == MODE_JPEG)
            save_jpg();
    }
    once=1;
}

int main (int argc, char ** argv)
{
    if (argc<2 || argc > 3)
    {  
        printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
        exit(0);
    }
    if(argc == 3)
    {
        mode = MODE_JPEG;
        filename = argv[2];
    }
    else if(argc == 2)
        mode = MODE_DISPLAY;
    
    glutInit(&argc,argv);
    loadScene(argv[2]);
    
    glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
    glutInitWindowPosition(0,0);
    glutInitWindowSize(WIDTH,HEIGHT);
    int window = glutCreateWindow("Ray Tracer");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    init();
    glutMainLoop();
}
