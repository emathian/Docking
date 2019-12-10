
//  EMILIE MATHIAN
// Compilation: g++  -o ARP_WARPQ1F ARP_WARPQ1F.cpp 
// Execution: ./RP_WARPQ1F test_q1.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stack>
#include <list>
#include <vector>
#include <queue>
#include <utility>
#include <iostream>

using namespace std;

#define MAX_ATOMS       10000
#define LINE_LENGTH     81
#define Pseudo_valence_angle_min 28
#define Pseudo_valence_angle_max 105
#define NUMATOMS 327

typedef struct {
    double  x, y, z;
} Point;


typedef struct Coord{
    int  i, j;
} Coord;

typedef struct Atom{
    int serial;
    char    atomName[5];
    char    altLoc[2];
    char    resName[4];
    char    chainID[2];
    int resSeq;
    char    iCode[2];
    Point   centre;
} Atom;


char* substring(const char* str, size_t begin, size_t len);

double dist_atom2(struct Atom *atom1, struct Atom *atom2);
/* Returns the Euclidian distance between two atoms. */
 
float pseudo_valence_angle(struct Atom *A, struct Atom *B, struct Atom *C);
/* Returns the angle according between three atoms.  */

/*
 * Declare an array to hold data read from the ATOM records of a PDB file.
 */

Atom    atom_c1[MAX_ATOMS+1];
Atom    atom_c2[MAX_ATOMS+1];
Coord   coord[MAX_ATOMS*5];


int read_data(char *filename, int n_file)
{
    FILE    *stream;
    char    line[LINE_LENGTH];
        char    s_serial[6];
        char    s_name[5];      /* Alpha-carbon is " CA "; calcium is "CA  " */
        char    s_altLoc[2];    /* Usually " " */
        char    s_resName[4];
        char    s_chainID[2];
        char    s_resSeq[5];
        char    s_iCode[2];     /* Usually " " */
        char    s_x[9];
        char    s_y[9];
        char    s_z[9];

        int     serial;
        int     resSeq;
        double  x;
        double  y;
        double  z;

    int i=0;

        if ( (stream = fopen(filename, "r")) == NULL ) {
                (void) fprintf(stderr, "Unable to open %s\n", filename);
                exit(0);
        }

        while ( fgets(line, LINE_LENGTH, stream) ) {
            if ( strncmp(line, "ATOM  ", 6) == 0 || strncmp(line, "HETATM", 6) == 0 ){
                
                 /* Split the line into its constituent fields.
                 * We are only interested in columns 1-54.
                 */

                strncpy(s_serial,  &line[6],  5); s_serial[5]  = '\0';
                strncpy(s_name,    &line[12], 4); s_name[4]    = '\0';
                strncpy(s_altLoc,  &line[16], 1); s_altLoc[1]  = '\0';
                strncpy(s_resName, &line[17], 3); s_resName[3] = '\0';
                strncpy(s_chainID, &line[21], 1); s_chainID[1] = '\0';
                strncpy(s_resSeq,  &line[22], 4); s_resSeq[4]  = '\0';
                strncpy(s_iCode,   &line[26], 1); s_iCode[1]   = '\0';
                strncpy(s_x,       &line[30], 8); s_x[8]       = '\0';
                strncpy(s_y,       &line[38], 8); s_y[8]       = '\0';
                strncpy(s_z,       &line[46], 8); s_z[8]       = '\0';

                /*
                 * Convert the numeric fields to integers or doubles.
                 * The library functions atoi() and atof() are
                 * described in the UNIX manual pages ('man atoi' and
                 * 'man atof').
                 */

                serial = atoi(s_serial);
                resSeq = atoi(s_resSeq);
                x      = atof(s_x);
                y      = atof(s_y);
                z      = atof(s_z);

                /*
                 * Copy values to the next element in the atom array.
                 */

                if ( ++i > MAX_ATOMS ) {
                    (void) fprintf(stderr, "Too many atoms read\n");
                    exit(0);
                }

                if (n_file == 1){
                    atom_c1[i].serial = serial;
                    strcpy(atom_c1[i].atomName, s_name);
                    strcpy(atom_c1[i].altLoc, s_altLoc);
                    strcpy(atom_c1[i].resName, s_resName);
                    strcpy(atom_c1[i].chainID, s_chainID);
                    atom_c1[i].resSeq = resSeq;
                    strcpy(atom_c1[i].iCode, s_iCode);
                    atom_c1[i].centre.x = x;
                    atom_c1[i].centre.y = y;
                    atom_c1[i].centre.z = z;
                }
                else{
                     atom_c2[i].serial = serial;
                    strcpy(atom_c2[i].atomName, s_name);
                    strcpy(atom_c2[i].altLoc, s_altLoc);
                    strcpy(atom_c2[i].resName, s_resName);
                    strcpy(atom_c2[i].chainID, s_chainID);
                    atom_c2[i].resSeq = resSeq;
                    strcpy(atom_c2[i].iCode, s_iCode);
                    atom_c2[i].centre.x = x;
                    atom_c2[i].centre.y = y;
                    atom_c2[i].centre.z = z;
                }     

      }
    }
    return i;
}

void write_pdb_atom(int serial,char *s_name, char *s_altLoc, char *s_resName, char *s_chainID,
      int  resSeq, char *s_iCode, Point centre)
{
    printf("ATOM  %5d %s%s%s %s%4d%s   %8.3f%8.3f%8.3f\n",
                                serial,
                                s_name,
                                s_altLoc,
                                s_resName,
                                s_chainID,
                                resSeq,
                                s_iCode,
                                centre.x,
                                centre.y,
                                centre.z);
}



int main(int argc, char ** argv)
{
    int numAtoms_f1;
    int numAtoms_f2;
    int i;
    int j;
    double threshold_min = 2.8;
    double threshold_max = 4.8;
    double d;

        if ( argc<=2 ) {
                (void) fprintf(stderr, "usage: atom_array file.pdb\n");
                exit(0);
        }

    numAtoms_f1 = read_data(argv[1], 1);

    numAtoms_f2 = read_data(argv[2], 2);

    for (i=1; i<=numAtoms_f2; ++i) {
        write_pdb_atom(
            atom_c2[i].serial,
            atom_c2[i].atomName,
            atom_c2[i].altLoc,
            atom_c2[i].resName,
            atom_c2[i].chainID,
            atom_c2[i].resSeq,
            atom_c2[i].iCode,
            atom_c2[i].centre);
    }


    return 0;
}


char* substring(const char* str, size_t begin, size_t len) 
{ 
  if (str == 0 || strlen(str) == 0 || strlen(str) < begin || strlen(str) < (begin+len)) 
    return 0; 

  return strndup(str + begin, len); 
} 


double dist_atom2(struct Atom *atom1, struct Atom *atom2) 
{ 
  /* Returns the Euclidian distance between two atoms. */

  double d  = sqrt(pow(atom1->centre.x - atom2->centre.x,2) + pow(atom1->centre.y - atom2->centre.y,2) + pow(atom1->centre.z - atom2->centre.z,2));
  return d;
} 


float pseudo_valence_angle(struct Atom *A, struct Atom *B, struct Atom *C)
{
  /* Returns the angle according between three atoms.  */

    // AB = B - A
    float AB[3];
    AB[0] = B->centre.x  - A->centre.x;
    AB[1] = B->centre.y  - A->centre.y;
    AB[2] = B->centre.z  - A->centre.z;
    // BC = C -B
    float BC[3];
    BC[0] = C->centre.x  - B->centre.x;
    BC[1] = C->centre.y  - B->centre.y;
    BC[2] = C->centre.z  - B->centre.z;


    // Dot product AB.AC
    float prod_scal = AB[0] * BC[0] + AB[1] * BC[1] + AB[2] * BC[2];

    // Norm AB
    float norm_AB = sqrt(pow(AB[0],2) + pow(AB[1],2) + pow(AB[2],2));
   
    // Norm AC
    float norm_BC = sqrt(pow(BC[0],2) + pow(BC[1],2) + pow(BC[2],2));

    // Angle
    double p =  prod_scal / (norm_AB * norm_BC);
    float theta = acos(p)*(180/M_PI); // Degree
    return theta;
}
