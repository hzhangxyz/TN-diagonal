#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <math.h>

typedef struct lattice{
    int n;
    int m;
    int point;
    int bonds_size;
    int (*bonds)[2];
    long matrix_size;
    float *matrix; // 1-H
    float *vector;
    float ans;
    int step;
} lattice;

int read_lattice(lattice* data, int argc, char** argv){
    data->n = atoi(argv[1]);
    data->m = atoi(argv[2]);
    data->point = data->n*data->m;
    data->matrix_size = 1;
    for(int i=0;i<data->point;i++)data->matrix_size*=2;
    data->bonds_size = data->n*(data->m-1) + data->m*(data->n-1);
    data->bonds = (int (*)[2])malloc(sizeof(int)*2*(data->n*(data->m-1)+data->m*(data->n-1)));
    data->vector = (float*)malloc(sizeof(float)*data->matrix_size);
    srand(time(NULL));
    for(int i=0;i<data->matrix_size;i++)data->vector[i]=((float)rand())/RAND_MAX;
    data->matrix = (float*)malloc(sizeof(float)*data->matrix_size*data->matrix_size);
    memset(data->matrix,0,sizeof(float)*data->matrix_size*data->matrix_size);
    for(long i=0;i<data->matrix_size;i++)data->matrix[(data->matrix_size+1)*i]=1;
#ifdef debug1
    for(int i=0;i<data->matrix_size;i++){
        for(int j=0;j<data->matrix_size;j++){
            printf("%.2f ",data->matrix[i*data->matrix_size+j]);
        }
        printf("\n");
    }
#endif
    return 0;
}

int generate_bond(lattice* data){
    for(int i=0;i<data->n-1;i++){
        for(int j=0;j<data->m;j++){
            data->bonds[i*data->m+j][0]=i*data->m+j;
            data->bonds[i*data->m+j][1]=(i+1)*data->m+j;
        }
    }
    for(int i=0;i<data->n;i++){
        for(int j=0;j<data->m-1;j++){
            data->bonds[data->m*(data->n-1)+i*(data->m-1)+j][0]=i*data->m+j;
            data->bonds[data->m*(data->n-1)+i*(data->m-1)+j][1]=i*data->m+j+1;
        }
    }
#ifdef debug0
    for(int i=0;i<(data->n*(data->m-1)+(data->m*(data->n-1)));i++)
        printf("%d\t%d\n",data->bonds[i][0],data->bonds[i][1]);
#endif
    return 0;
}

int _set_matrix(lattice* data, int a, int b, int now, long offset_x, long offset_y, float value, int flag){
    if(now==data->point){
        data->matrix[offset_x*data->matrix_size+offset_y] -= value;
        return 0;
    }
    if(now==b){
        switch(flag){
            case(0):
                _set_matrix(data, a, b, now+1, offset_x*2  , offset_y*2,   +0.25, 0);
                _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+1, -0.25, 0);
                break;
            case(1):
                _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2,   0.5,   0);
                break;
            case(2):
                _set_matrix(data, a, b, now+1, offset_x*2,   offset_y*2+1, 0.5,   0);
                break;
            case(3):
                _set_matrix(data, a, b, now+1, offset_x*2  , offset_y*2,   -0.25, 0);
                _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+1, +0.25, 0);
                break;
        }
    }else if(now==a){
        _set_matrix(data, a, b, now+1, offset_x*2  , offset_y*2,   value, 0);
        _set_matrix(data, a, b, now+1, offset_x*2  , offset_y*2+1, value, 1);
        _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2,   value, 2);
        _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+1, value, 3);
    }else{
        _set_matrix(data, a, b, now+1, offset_x*2,   offset_y*2,   value, flag);
        _set_matrix(data, a, b, now+1, offset_x*2+1, offset_y*2+1, value, flag);
    }
    return 0;
}

int set_matrix(lattice* data){
    for(int bond=0;bond<data->bonds_size;bond++){
        _set_matrix(data,data->bonds[bond][0],data->bonds[bond][1],0,0,0,0,0);
    }
#ifdef debug2
    for(int i=0;i<data->matrix_size;i++){
        for(int j=0;j<data->matrix_size;j++){
            if(data->matrix[i*data->matrix_size+j]==0.)
                printf("      ");
            else
                printf("%5.2f ",data->matrix[i*data->matrix_size+j]);
        }
        printf("\n");
    }
#endif
    return 0;
}

int _update_vector(lattice* data){
    float* temp = (float*)malloc(sizeof(float)*data->matrix_size);
    for(int i=0;i<data->matrix_size;i++){
        temp[i]=data->vector[i];
        data->vector[i]=0;
    }
    for(int i=0;i<data->matrix_size;i++){
        for(int j=0;j<data->matrix_size;j++){
            data->vector[i] += data->matrix[data->matrix_size*i+j]*temp[j];
        }
    }
    free(temp);
    return 0;
}

int find_max(lattice *data){
#ifdef debug3
    for(int i=0;i<data->matrix_size;i++)printf("%5.2f ",data->vector[i]);
    printf("\n");
#endif
    data->ans = 0;
    int flag = 1;
    for(int t=0;(t<10)||(flag);t++){
        float temp=0;
        for(int i=0;i<data->matrix_size;i++){
            temp+=data->vector[i]*data->vector[i];
        }
        temp = sqrt(temp);
#define MIN 0.00000001
        if((data->ans-temp<MIN)&&(temp-data->ans<MIN)){
            flag=0;
            data->step = t;
        }
        data->ans = temp;
        for(int i=0;i<data->matrix_size;i++){
            data->vector[i] /= data->ans;
        }
        _update_vector(data);
    }
    return 0;
}

int output(lattice *data){
    printf("Size:\t%d\t%d\nStep:\t%d\nEnergy:\t%.6f\n\n",data->n,data->m,data->step,(1-data->ans)/(data->n*data->m));
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

