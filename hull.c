#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#define BUFSIZE 100

// Structure du Vecteur
struct vec {
    double x;
    double y;
};


// Fonction pour calculer le produit scalaire de deux vecteurs
// Formule : (x1 * x2) + (y1 * y2)

double dot(const struct vec *v1, const struct vec *v2) {
    return v1->x * v2->x + v1->y * v2->y;
}

// Fonction pour calculer la norme d'un vecteur
// Formule : (x2 − x1)(y3 − y1) − (y2 − y1)(x3 − x1)
double cross(const struct vec *p1,const struct vec *p2, const struct vec *p3){
    return (p2->x - p1->x)*(p3->y - p1->y) - (p2->y - p1->y)*(p3->x - p1->x);
}

// Fonction qui regarde si la suite de point constitue un tournant a gauche
// Formule : si le produit vectoriel est positif alors P3 est a gauche de P1P2
// et negatif si a droite
bool is_left_turn(const struct vec *p1,const struct vec *p2, const struct vec *p3){
    return cross(p1,p2,p3) > 0;
}

struct vecset {
    struct vec *data;
    size_t size;
    size_t capacity;
};

void vecset_create(struct vecset *self){
    self->capacity = 1;
    self->size = 0;
    self->data = malloc(self->capacity * sizeof(struct vec));
}

void vecset_destroy(struct vecset *self){
    free(self->data);
    self->size = 1;
    self->capacity = 1;
    self->data = NULL;
}

void vecset_add(struct vecset *self, struct vec p){
    if(self->size >= self->capacity) {
        self->capacity = self->capacity * 2;
        struct vec *data2 = malloc(self->capacity * sizeof(struct vec));
        data2 = memcpy(data2, self->data, self->size * sizeof(struct vec));
        free(self->data);
        self->data = data2;
    }
    self->size = self->size + 1;
    self->data[self->size - 1] = p;
}

//On considère une fonction de comparaison de points avec un contexte qui
//renvoie un entier strictement négatif si p1 est «plus petit» que p2, un entier
//strictement positif si p1 est «plus grand» que p2 et 0 si p1 est «égal» à p2.
typedef int (*comp_func_t)(const struct vec *p1,const struct vec *p2, const void *ctx);

int cmp1(const struct vec *p1,const struct vec *p2, const void *ctx){
    if(p1->y < p2->y){
        return -1;
    }
    if(p1->y > p2->y){
        return 1;
    }
    else {
        return 0;
    }
}

typedef int (*comp_func_t)(const struct vec *p1, const struct vec *p2, const void *ctx);

const struct vec *vecset_max(const struct vecset *self, comp_func_t func, const void *ctx){
    int max_idx = 0;
    for(int i = 1; i < self->size; i++) {
        if(func(&self->data[max_idx], &self->data[i], ctx) < 0) {
            max_idx = i;
        }
    }
    return &self->data[max_idx];
}


//ode d’une fonction qui renvoie le minimum d’un
//ensemble de points suivant une fonction de comparaison donnée
const struct vec *vecset_min(const struct vecset *self,comp_func_t func, const void *ctx){
    int min_idx = 0;
    for(int i = 1; i < self->size; i++) {
        if(func(&self->data[min_idx], &self->data[i], ctx) > 0) {
            min_idx = i;
        }
    }
    return &self->data[min_idx];
}


// Donner le code d’une fonction qui trie l’ensemble de points
//suivant la fonction de comparaison donnée
//tri par selection mais on peut faire mieux ( a faire plus tard)
void vecset_sort(struct vecset *self, comp_func_t func,const void *ctx){
    for(int i = 0; i < self->size - 1; i ++){
        for(int j = i + 1; j < self->size; j ++){
            if(func(&self->data[i], &self->data[j], ctx) > 0){
                struct vec temp = self->data[i];
                self->data[i] = self->data[j];
                self->data[j] = temp;
            }
        }
    }

}

// Fonction qui empile un élément.
void vecset_push(struct vecset *self, struct vec p){
    vecset_add(self, p);
}

//fonction qui dépile un élément.
void vecset_pop(struct vecset *self){
    self->size = self->size - 1;
    if(self->size <= self->capacity/4) {
        self->data = realloc(self->data, self->capacity / 2 * sizeof(struct vec));
        self->capacity = self->capacity / 2;
    }

}

//code d’une fonction qui renvoie le premier élément
//de la pile.
const struct vec *vecset_top(const struct vecset *self){
    return &self->data[self->size - 1];
}

//code d’une fonction qui renvoie le second élément
//de la pile.
const struct vec *vecset_second(const struct vecset *self){
    return &self->data[self->size - 2];
}


void jarvis_march(const struct vecset *in, struct vecset *out){

}



int main() {
    //setbuf(stdout, NULL); // avoid buffering in the output

    char buffer[BUFSIZE];
    fgets(buffer, BUFSIZE, stdin);

    size_t count = strtol(buffer, NULL, 10);

    struct vecset *self = malloc(sizeof(struct vecset));
    vecset_create(self);
    printf("%d, %d\n", self->size, self->capacity);


    for (size_t i = 0; i < count; ++i) {
        struct vec p;

        fgets(buffer, BUFSIZE, stdin);

        char *endptr = buffer;
        p.x = strtod(endptr, &endptr);
        p.y = strtod(endptr, &endptr);

        // then do something with p and test function
        printf("%f %f\n", p.x, p.y);
        vecset_add(self, p);


    }


    //test des fonction avec struct vec p
    printf("\nMaximum\n");
    printf("%f %f\n", vecset_max(self, cmp1, NULL)->x, vecset_max(self, cmp1, NULL)->y);

    printf("\nMinimum\n");
    printf("%f %f\n", vecset_min(self, cmp1, NULL)->x, vecset_min(self, cmp1, NULL)->y);

    printf("\nTri\n");
    vecset_sort(self, cmp1, NULL);
    for(int i = 0; i < self->size; i++){
        printf("%f %f\n", self->data[i].x, self->data[i].y);
    }

    printf("\nEmpile\n");
    struct vec p;
    p.x = 1;
    p.y = 2;
    vecset_push(self, p);

    printf("\nDepile\n");
    vecset_pop(self);

    printf("\nTop\n");
    printf("%f %f\n", vecset_top(self)->x, vecset_top(self)->y);

    printf("\nSecond\n");
    printf("%f %f\n", vecset_second(self)->x, vecset_second(self)->y);

    printf("\ndot\n");
    printf("%f\n", dot(&self->data[0], &self->data[1]));

    printf("\ncross\n");
    printf("%f\n", cross(&self->data[0], &self->data[1], &self->data[2]));

    printf("\nis_left_turn\n");
    printf("%d\n", is_left_turn(&self->data[0], &self->data[1], &self->data[2]));







    vecset_destroy(self);
    return 0;
}

//Exemple execution
//./hull < input.txt > output.txt 2>hull.log
//chmod +x ./hull-generator