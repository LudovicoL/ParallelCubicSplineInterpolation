/* Generate an n-element double vector and write it to a file */

#include <stdio.h>
#include <stdlib.h>

#define MAX 2
#define MIN 1


int main (int argc, char * argv[]) {
    int i, sign;
    int n;
    FILE *foutptr;
    double *array;
    if(argc != 3) {
        printf("Error! To execute:\n./gen number_of_elements output_filename\n");
        return -1;
    }

    printf ("argv[0] is '%s'\n", argv[0]);
    printf ("argv[1] is '%s'\n", argv[1]);
    printf ("argv[2] is '%s'\n", argv[2]);
    n = atoi (argv[1]);
    array = (double *) malloc (n * sizeof(double));
    foutptr = fopen (argv[2], "w");

    fprintf (foutptr, "%d\n", n);       // Save value n

    array[0] = 1;
    for (i = 1; i <= n; i++) {          // Generate and save elements x
        array[i] = array[i-1] + (double) (rand() % (MAX - MIN + 1) + MIN);
        fprintf (foutptr, "%lf ", array[i]);
    }

    fprintf (foutptr, "\n");

    for (i = 0; i < n; i++) {           // Generate and save elements y
        do {
            sign = (rand() % (2));
		    if(sign == 0) {sign -= 1;}
            array[i] = ((double)(i) / (double)(MAX+MIN)) + sign * (rand() % (MAX - MIN + 1) + MIN);
        } while(array[i] < 0);
        fprintf (foutptr, "%lf ", array[i]);
    }
    
    fclose (foutptr);
    return 0;
}
