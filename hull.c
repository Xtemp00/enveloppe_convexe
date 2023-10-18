#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#define BUFSIZE 100

// Structure du Vecteur
struct vec {
    double x;
    double y;
};


// Fonction pour calculer la distance entre deux points
// Formule : (x1 * x2) + (y1 * y2)

double dot(const struct vec *v1, const struct vec *v2) {
    return v1->x * v2->x + v1->y * v2->y;
}

// Fonction pour calculer la norme d'un vecteur
// Formule : (x2 − x1)(y3 − y1) − (y2 − y1)(x3 − x1)

double cross(const struct vec *p1,const struct vec *p2, const struct vec *p3){
    return (p2->x - p1->x)*(p3->y - p1->y) - (p2->y - p1->y)*(p3->x - p1->x);
}
// il me casse les coilles a tout le temps parler
bool is_left_turn(const struct vec *p1,const struct vec *p2, const struct vec *p3){

}

struct vecset {
    struct vec *data;
    size_t size;
    size_t capacity;
};

void vecset_create(struct vecset *self){
    self->capacity = 0;
    self->size = 0;
    self->data = NULL;
}

void vecset_destroy(struct vecset *self){
    free(self->data);
    self->size = 0;
    self->capacity = 0;
    self->data = NULL;
}

void vecset_add(struct vecset *self, struct vec p){
    self->size = self->size + 1;
    if(self->size >= self->capacity) {
        self->data = realloc(self->data, self->capacity * 2 * sizeof(struct vec));
        self->capacity = self->capacity * 2;
    }
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
    struct vec max = self->data[0];
    for(int i = 0; i < self->size - 1; i ++){
        if(func(&self->data[i], &self->data[i+1], ctx) > 0){
            max = self->data[i];
        }
    }
}

//ode d’une fonction qui renvoie le minimum d’un
//ensemble de points suivant une fonction de comparaison donnée
const struct vec *vecset_min(const struct vecset *self,comp_func_t func, const void *ctx){
    struct vec min = self->data[0];
    for(int i = 0; i < self->size - 1; i ++){
        if(func(&self->data[i], &self->data[i+1], ctx) < 0){
            min = self->data[i];
        }
    }

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

void vecset_push(struct vecset *self, struct vec p){
    
}

void vecset_pop(struct vecset *self){

}


int main() {
    //setbuf(stdout, NULL); // avoid buffering in the output

    char buffer[BUFSIZE];
    fgets(buffer, BUFSIZE, stdin);

    size_t count = strtol(buffer, NULL, 10);

    for (size_t i = 0; i < count; ++i) {
        struct vec p;

        fgets(buffer, BUFSIZE, stdin);

        char *endptr = buffer;
        p.x = strtod(endptr, &endptr);
        p.y = strtod(endptr, &endptr);

    // then do something with p
    }
    return 0;
}

//Exemple execution
//./hull < input.txt > output.txt 2>hull.log