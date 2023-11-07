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
    self->data[self->size] = p;
    self->size = self->size + 1;
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

//farthest point from line (XY )

void farthest_point(const struct vecset *in, struct vec *out, const struct vec *X, const struct vec *Y){
    double max = 0;
    for(int i = 0; i < in->size; i++){
        double distance = fabs(cross(X, Y, &in->data[i]));
        if(distance > max){
            max = distance;
            *out = in->data[i];
        }
    }
}

//function JarvisMarch(S)
//R←∅
//F ← leftmost point in S
//C←F
//repeat
//R ← R ∪ {C}
//N ← a point in S
//for all I ∈ S do
//if (C, I, N ) is a left turn then
//N ←I
//end if
//end for
//C←N
//until F = C
//return R
//end function
//Marche de Jarvis
void jarvis_march(const struct vecset *in, struct vecset *out) {
    // On commence par trouver le point le plus a gauche
    const struct vec *F = vecset_min(in, cmp1, NULL);
    // On l'ajoute a l'enveloppe convexe
    vecset_add(out, *F);
    // On initialise le point courant
    const struct vec *C = F;
    // On initialise le point suivant
    const struct vec *N = NULL;
    // On initialise le point suivant
    do {
        N = &in->data[0];
        // On parcourt tous les points
        for(int i = 1; i < in->size; i++){
            // Si le point est a gauche du point courant
            if(is_left_turn(C, N, &in->data[i])){
                // On le met dans le point suivant
                N = &in->data[i];
            }
        }
        // On ajoute le point suivant a l'enveloppe convexe
        vecset_add(out, *N);
        // On met le point suivant dans le point courant
        C = N;
    } while (C != F);
    // il manque le dernier point
    vecset_pop(out);
}



// Algorithm 2 Parcours de Graham
// function (S)
// B ← lowest point in S
// Sort(S)
// F ← first point in S
// Push(R, B)
// Push(R, F )
// for all I ∈ S \ {B, F } do
// T ← top point of R
// S ← second point of R
// while |R| ≥ 2 and (S, T, I) is a left turn do
// Pop(R)
// end while
// Push(R, I)
// end for
// return R
// end function
/*
// Fonction de comparaison utilisant atan2
int cmp_angle(const void *a, const void *b) {
    const struct vec *p1 = (const struct vec *)a;
    const struct vec *p2 = (const struct vec *)b;
    double angle1 = atan2(p1->y - ref_point.y, p1->x - ref_point.x);
    double angle2 = atan2(p2->y - ref_point.y, p2->x - ref_point.x);

    if (angle1 < angle2) return -1;
    if (angle1 > angle2) return 1;

    // Si les angles sont égaux, le point le plus proche de ref_point devrait venir en premier
    double distance1 = (p1->x - ref_point.x) * (p1->x - ref_point.x) + (p1->y - ref_point.y) * (p1->y - ref_point.y);
    double distance2 = (p2->x - ref_point.x) * (p2->x - ref_point.x) + (p2->y - ref_point.y) * (p2->y - ref_point.y);

    return (distance1 > distance2) ? -1 : (distance1 < distance2);
}

// Graham Scan
void graham_scan(const struct vecset *in, struct vecset *out) {
    // Trouver le point avec l'ordonnée la plus basse et le mettre en tant que ref_point.
    ref_point = *vecset_min(in, (comp_func_t)cmp1, NULL);

    // Trier les points par angle polaire par rapport à 'ref_point'.
    qsort(in->data, in->size, sizeof(struct vec), cmp_angle);

    // Initialiser la pile 'out' avec le point le plus bas.
    vecset_push(out, ref_point);

    // Parcourir les points triés et construire l'enveloppe convexe.
    for (size_t i = 1; i < in->size; ++i) {
        // Gardez un œil sur le sommet actuel de la pile.
        while (out->size >= 2 && !is_left_turn(vecset_second(out), vecset_top(out), &in->data[i])) {
            // Si ce n'est pas un tournant à gauche, retirez le point du sommet.
            vecset_pop(out);
        }

        // Ajoutez le prochain point à la pile.
        vecset_push(out, in->data[i]);
    }
}

*/


//function FindHull(S, X, Y )
//if S = ∅ then
//return ∅
//end if
//M ← farthest point from line (XY )
//S1 ← ∅
//S2 ← ∅
//for all I ∈ S \ {M } do
//−−→
//if I in on the left of XM then
//S1 ← S1 ∪ {I}
//end if
//−−→
//if I in on the left of M Y then
//S2 ← S2 ∪ {I}
//end if
//end for
//R1 ← FindHull(S1 , X, M )
//R2 ← FindHull(S2 , M, Y )
//return R1 ∪ {M } ∪ R2
//end function

void FindHull(const struct vecset *in, struct vecset *out, const struct vec *X, const struct vec *Y){
    if(in->size == 0){
        return;
    }
    const struct vec *M = malloc(sizeof(struct vec));
    farthest_point(in, M, X, Y);

    struct vecset *S1 = malloc(sizeof(struct vecset));
    vecset_create(S1);


    struct vecset *S2 = malloc(sizeof(struct vecset));
    vecset_create(S2);

    for(int i = 0; i < in->size; i++){
        if((is_left_turn(X, M, &in->data[i])) && (&in->data[i] != M)){
            vecset_add(S1, in->data[i]);
        }
        if((is_left_turn(M, Y, &in->data[i])) && (&in->data[i] != M)){
            vecset_add(S2, in->data[i]);
        }
    }

    struct vecset *R1 = malloc(sizeof(struct vecset));
    vecset_create(R1);

    struct vecset *R2 = malloc(sizeof(struct vecset));
    vecset_create(R2);

    FindHull(S1, R1, X, M);
    FindHull(S2, R2, M, Y);

    for(int i = 0; i < R1->size; i++){
        vecset_add(out, R1->data[i]);
    }
    vecset_add(out, *M);
    for(int i = 0; i < R2->size; i++){
        vecset_add(out, R2->data[i]);
    }

    free(S1);
    free(S2);
    free(R1);
    free(R2);

}

//function QuickHull(S)
//A ← leftmost point in S
//B ← rightmost point in S
//S1 ← ∅
//S2 ← ∅
//for all I ∈ S \ {A, B} do
//−−→
//if I is on the left of AB then
//S1 ← S1 ∪ {I}
//else
//S2 ← S2 ∪ {I}
//end if
//end for
//R1 ← FindHull(S1 , A, B)
//R2 ← FindHull(S2 , B, A)
//return {A} ∪ R1 ∪ {B} ∪ R2
//end function

void quickhull(const struct vecset *in, struct vecset *out){
    // On commence par trouver le point le plus a gauche
    const struct vec *A = vecset_min(in, cmp1, NULL);
    // On trouve le point le plus a droite
    const struct vec *B = vecset_max(in, cmp1, NULL);
    struct vecset *S1 = malloc(sizeof(struct vecset));
    vecset_create(S1);

    struct vecset *S2 = malloc(sizeof(struct vecset));
    vecset_create(S2);

    for(int i = 0; i < in->size; i++){
        if(is_left_turn(A, B, &in->data[i]) && (&in->data[i] != A) && (&in->data[i] != B)){
            vecset_add(S1, in->data[i]);
        }
        else{
            vecset_add(S2, in->data[i] );
        }
    }

    struct vecset *R1 = malloc(sizeof(struct vecset));
    vecset_create(R1);

    struct vecset *R2 = malloc(sizeof(struct vecset));
    vecset_create(R2);

    FindHull(S1, R1, A, B);
    FindHull(S2, R2, B, A);
    vecset_add(out, *A);
    for(int i = 0; i < R1->size; i++){
        vecset_add(out, R1->data[i]);
    }
    vecset_add(out, *B);
    for(int i = 0; i < R2->size; i++){
        vecset_add(out, R2->data[i]);
    }

    free(S1);
    free(S2);
    free(R1);
    free(R2);

}




int main() {
    //setbuf(stdout, NULL); // avoid buffering in the output

    char buffer[BUFSIZE];
    fgets(buffer, BUFSIZE, stdin);

    size_t count = strtol(buffer, NULL, 10);

    struct vecset *self = malloc(sizeof(struct vecset));
    vecset_create(self);
    //printf("%d, %d\n", self->size, self->capacity);


    for (size_t i = 0; i < count; ++i) {
        struct vec p;

        fgets(buffer, BUFSIZE, stdin);

        char *endptr = buffer;
        p.x = strtod(endptr, &endptr);
        p.y = strtod(endptr, &endptr);

        // then do something with p and test function
        //printf("%f %f\n", p.x, p.y);
        vecset_add(self, p);


    }

    /*
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


    //test de jarvis march
    struct vecset *out = malloc(sizeof(struct vecset));
    vecset_create(out);
    jarvis_march(self, out);
    // On Affiche la taille en début de fichier
    printf("%li\n", out->size);
    // On print les points de l'enveloppe convexe
    for(int i = 0; i < out->size; i++){
        printf("%f %f\n", out->data[i].x, out->data[i].y);
    }


    //test de graham
    struct vecset *out = malloc(sizeof(struct vecset));
    vecset_create(out);
    graham_scan(self, out);
    // On Affiche la taille en début de fichier
    printf("%li\n", out->size);
    // On print les points de l'enveloppe convexe
    for(int i = 0; i < out->size; i++){
        printf("%f %f\n", out->data[i].x, out->data[i].y);
    }
    */

    //Quick Hull

    //test de jarvis march
    struct vecset *out = malloc(sizeof(struct vecset));
    vecset_create(out);
    quickhull(self, out);
    // On Affiche la taille en début de fichier
    printf("%li\n", out->size);
    // On print les points de l'enveloppe convexe
    for(int i = 0; i < out->size; i++){
        printf("%f %f\n", out->data[i].x, out->data[i].y);
    }


    vecset_destroy(self);
    return 0;
}

//Exemple execution
//./hull < input.txt > output.txt 2>hull.log
//chmod +x ./hull-generator
