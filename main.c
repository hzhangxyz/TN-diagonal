#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

typedef struct matrix_element{
    long                   x;
    long                   y;
    float                  value;
    struct matrix_element* left;
    struct matrix_element* right;
} matrix_element;

typedef struct lattice{
    int             n;
    int             m;
    int             point;
    int             bonds_size;
    int             (*bonds)[2];
    long            matrix_size;
    matrix_element* matrix;
    float           *vector;
    float           ans;
    int             step;
    double          time_start;
    double          time_end;
} lattice;

inline int init_element(matrix_element** matrix, long x, long y, float value){
    *matrix          = (matrix_element*)malloc(sizeof(matrix_element));
    (*matrix)->x     = x;
    (*matrix)->y     = y;
    (*matrix)->value = value;
    (*matrix)->left  = NULL;
    (*matrix)->right = NULL;
    return 0;
}

int add_element(matrix_element* matrix, long x, long y, float value){
    if((matrix->x==x)&&(matrix->y==y)){
        matrix->value += value;
    }else if((matrix->x<x)||((matrix->x==x)&&(matrix->y<y))){
        if(matrix->left){
            add_element(matrix->left,x,y,value);
        }else{
            init_element(&matrix->left,x,y,value);
        }
    }else{
        if(matrix->right){
            add_element(matrix->right,x,y,value);
        }else{
            init_element(&matrix->right,x,y,value);
        }
    }
    return 0;
}

int traversal(matrix_element* matrix, float* src, float* dst){
    dst[matrix->y] += matrix->value*src[matrix->x];
    if(matrix->left){
        traversal(matrix->left,src,dst);
    }
    if(matrix->right){
        traversal(matrix->right,src,dst);
    }
    return 0;
}

inline int read_lattice(lattice* data, int argc, char** argv){
    data->time_start  = (double)clock()/CLOCKS_PER_SEC;
    data->n           = atoi(argv[1]);
    data->m           = atoi(argv[2]);
    data->point       = data->n*data->m;
    data->matrix_size = 1;
    for(int i=0;i<data->point;i++){
        data->matrix_size *= 2;
    }
    data->bonds_size  = data->n*(data->m-1) + data->m*(data->n-1);
    data->bonds       = (int (*)[2])malloc(sizeof(int)*2*(data->n*(data->m-1)+data->m*(data->n-1)));
    data->vector      = (float*)malloc(sizeof(float)*data->matrix_size);
    srand(time(NULL));
    for(long i=0;i<data->matrix_size;i++){
        data->vector[i]   = ((float)rand())/RAND_MAX;
    }
    init_element(&data->matrix,0,0,0);
    for(long i=0;i<data->matrix_size;i++){
        add_element(data->matrix,i,i,1);
    }
    return 0;
}

inline int generate_bond(lattice* data){
    for(int i=0;i<data->n-1;i++){
        for(int j=0;j<data->m;j++){
            data->bonds[i*data->m+j][0] = i*data->m+j;
            data->bonds[i*data->m+j][1] = (i+1)*data->m+j;
        }
    }
    for(int i=0;i<data->n;i++){
        for(int j=0;j<data->m-1;j++){
            data->bonds[data->m*(data->n-1)+i*(data->m-1)+j][0] = i*data->m+j;
            data->bonds[data->m*(data->n-1)+i*(data->m-1)+j][1] = i*data->m+j+1;
        }
    }
    return 0;
}

int _set_matrix(lattice* data, int a, int b, int now, long offset_x, long offset_y, float value, int flag){
    if(now==data->point){
        add_element(data->matrix,offset_x,offset_y,-value);
        return 0;
    }
    if(now==b){
        switch(flag){
            case(0):
                _set_matrix(data, a, b, now+1, offset_x*2+0, offset_y*2+0, +0.25, 0);
                _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+1, -0.25, 0);
                break;
            case(1):
                _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+0, +0.50, 0);
                break;
            case(2):
                _set_matrix(data, a, b, now+1, offset_x*2+0, offset_y*2+1, +0.50, 0);
                break;
            case(3):
                _set_matrix(data, a, b, now+1, offset_x*2+0, offset_y*2+0, -0.25, 0);
                _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+1, +0.25, 0);
                break;
        }
    }else if(now==a){
        _set_matrix(data, a, b, now+1, offset_x*2+0, offset_y*2+0, value, 0);
        _set_matrix(data, a, b, now+1, offset_x*2+0, offset_y*2+1, value, 1);
        _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+0, value, 2);
        _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+1, value, 3);
    }else{
        _set_matrix(data, a, b, now+1, offset_x*2+0, offset_y*2+0, value, flag);
        _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+1, value, flag);
    }
    return 0;
}

inline int set_matrix(lattice* data){
    for(int bond=0;bond<data->bonds_size;bond++){
        _set_matrix(data,data->bonds[bond][0],data->bonds[bond][1],0,0,0,0,0);
    }
    return 0;
}

inline int _update_vector(lattice* data){
    float* temp = (float*)malloc(sizeof(float)*data->matrix_size);
    for(long i=0;i<data->matrix_size;i++){
        temp[i]         = data->vector[i];
        data->vector[i] = 0;
    }
    traversal(data->matrix,temp,data->vector);
    free(temp);
    return 0;
}

#define MIN 1e-8
inline int find_max(lattice *data){
    data->ans = 0;
    for(int t=0;;t++){
        float temp=0;
        for(long i=0;i<data->matrix_size;i++){
            temp            += data->vector[i]*data->vector[i];
        }
        temp = sqrt(temp);
        if((data->ans-temp<MIN)&&(temp-data->ans<MIN)&&(t>10)){
            data->step      = t;
            return 0;
        }
        data->ans = temp;
        for(long i=0;i<data->matrix_size;i++){
            data->vector[i] /= data->ans;
        }
        _update_vector(data);
    }
}

inline int output(lattice *data){
    data->time_end = (double)clock() / CLOCKS_PER_SEC;
    printf(
            "Size:\t%d\t%d\nTime:\t%.6f\nStep:\t%d\nEnergy:\t%.6f\nValue:\n",
            data->n,
            data->m,
            data->time_end-data->time_start,
            data->step,
            (1-data->ans)/(data->n*data->m)
    );
    for(long i=0;i<data->matrix_size;i++){
        printf("%.2f ",data->vector[i]);
    }
    printf("\n\n");
    free(data->bonds);
    free(data->vector);
    free(data->matrix);
    return 0;
}

int main(int argc, char** argv){
    lattice data;
    read_lattice(&data,argc,argv);
    generate_bond(&data);
    set_matrix(&data);
    find_max(&data);
    output(&data);
}

