#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

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

bool is_left_turn(const struct vec *p1,const struct vec *p2, const struct vec *p3){

}
