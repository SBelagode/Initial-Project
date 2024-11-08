#include <iostream>
//import library for computing interval math
#include "mpfi.h"
#include "mpfi_io.h"
#include <queue>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>



/*PLEASE READ: 
Sources: https://math.libretexts.org/Bookshelves/Calculus/CLP-4_Vector_Calculus_(Feldman_Rechnitzer_and_Yeager)/06%3A_Appendices/6.01%3A_A_Appendices/6.1.08%3A_A.8_Conic_Sections_and_Quadric_Surfaces
https://www.staff.city.ac.uk/o.castro-alvaredo/teaching/surfaces.pdf 

Stucture of Program: 
- Defined Functions: 
    --intervalClass: 
        ---Uses mpfi and mpfr libraries to define l x w x h boxes that I use to compute whether both g and f are defined. This 
           makes it easier to use within the queue structure. These boxes also define methods to print the intervals, and constructors
           to initialize with either explicit double values or copy from a previously defined interval box (this latter case is especially useful 
           as the precision of mpfr extends beyond double precision). There are also methods to define whether or not the max diameter
           of any of the intervals is less than or equal to a particular precision value. Until I gain a better understanding of the mpfr 
           library, the precision will only be that of a double
    -- f (sphere) and g (elliptic paraboloid)
        ---These are simply just the functions f = x^2 + y^2 + z^2 -1 and g = x^2 + y^2 - z - 1 with arguments being the x,y,and 
           z intervals
    -- hasZero
        ---This function just computes f and g for a given x,y,z interval box and checks if both resulting intervals contain 0


    -- hyperboloid of one sheet: x^2/a^2 + y^2/b^2 - z^2/c^2 = 1

    -- hyperboloid of two sheets: x^2/a^2 + y^2/b^2 - z^2/c^2 = -1

    -- hyperplane: aX + bY + cZ + d = 0

    -- ellipsoid: x^2/a^2 + y^2/b^2 + z^2/c^2 = 1

    -- elliptic cone: x^2 / a^2 + y^2 / b^2 = z^2 /c^2

    -- hyperbolic paraboloid: x^2/a^2 - y^2/b^2 = z/c

    -- elliptic paraboloid: x^2/a^2 + y^2/b^2 - z = 0 

    The following are boolean valued functions returning true if the equation evaluated at the given intervals returns 0 or not (these are mainly conic sections):
    -- elliptic cylinder: x^2/a + y^2/b = 1

    -- parabolic cylinder: y = ax^2

    -- hyperbolic cylinder: x^2/ a - y^2/ b = 1
    
- main:
    --Define a maxError double value. This is essentially the desired width threshold that one would like intervals 
      that define the boxes to be within 
    -- Define initial box to compute within 
    -- define queue and vector (vector holds the boxes that are smaller than the maxError and where f and g contain 0), queue initial
       box
    -- while queue is not empty, process and pop first interval box, check if f and g contain 0 inside this box. If no, process next 
       item in queue. If yes, then check if box is smaller than maxError threshold. if so, add the box to the answer vector. if 
       not, then split box into 8 smaller identical boxes and queue all 8. 
    -- outside while loop print all the boxes in the answer vector and clear memory used. 
- testing code is commented outside of main
- manually assign error, than output of program is list of boxes

- can verify graphs with https://www.geogebra.org/3d 

*/




/*
Creating intervalClass that basically serves as a blueprint to create interval objects. So each object would represent
a box interval in which we will evaluate f and g. This is done so that we can work with these interval structs in a queue, 
as when I was trying to queue these intervals by themselves, I was not able to. 
*/

class intervalClass {
    public:
    mpfi_t  x;
    mpfi_t  y;
    mpfi_t  z;
    mpfr_t maxDiameterX;
    mpfr_t maxDiameterY;
    mpfr_t maxDiameterZ;

    mpfr_t intermediateMax1;
    mpfr_t intermediateMax2;

    //default constructor
    intervalClass(){

    }

    // initialize variables and set them equal to given bounds
    intervalClass(double x1, double x2, double y1, double y2, double z1, double z2) {
        mpfi_init((mpfi_ptr)x);
        mpfi_interv_d((mpfi_ptr)x, x1, x2);
        mpfi_init((mpfi_ptr)y);
        mpfi_interv_d((mpfi_ptr)y, y1, y2);
        mpfi_init((mpfi_ptr)z);
        mpfi_interv_d((mpfi_ptr)z, z1, z2);

        mpfr_init(maxDiameterX);
        mpfr_init(maxDiameterY);
        mpfr_init(maxDiameterZ);

        mpfr_init(intermediateMax1);
        mpfr_init(intermediateMax2);
    }

    //initialize from existing values
    intervalClass(mpfi_t xcopy, mpfi_t ycopy, mpfi_t zcopy) {
        mpfi_init((mpfi_ptr)x);
        mpfi_set((mpfi_ptr)x, xcopy);
        mpfi_init((mpfi_ptr)y);
        mpfi_set((mpfi_ptr)y, ycopy);
        mpfi_init((mpfi_ptr)z);
        mpfi_set((mpfi_ptr)z, zcopy);

        mpfr_init(maxDiameterX);
        mpfr_init(maxDiameterY);
        mpfr_init(maxDiameterZ);

        mpfr_init(intermediateMax1);
        mpfr_init(intermediateMax2);
    }


    

    //print method

    void p(FILE *strm){
        fprintf(strm, "X values: \n");
        mpfi_out_str(strm, 10, 0, (mpfi_ptr) x);
        fprintf(strm, "\n");
        fprintf(strm, "Y values: \n");
        mpfi_out_str(strm, 10, 0, (mpfi_ptr) y);
        fprintf(strm, "\n");
        fprintf(strm, "Z values: \n");
        mpfi_out_str(strm, 10, 0, (mpfi_ptr) z);
        fprintf(strm, "\n");
        fprintf(strm, "\n");
    }

    void p(){
        std::cout << "X values: " << "\n";
        mpfi_out_str(stdout, 10, 0, (mpfi_ptr) x);
        std::cout << "\n";
        std::cout << "Y values: " << "\n";
        mpfi_out_str(stdout, 10, 0, (mpfi_ptr) y);
        std::cout << "\n";
        std::cout << "Z values: " << "\n";
        mpfi_out_str(stdout, 10, 0, (mpfi_ptr) z);
        std::cout << "\n";
        std::cout << "\n";
    }

    //gets max diameter and then checks if it is less than or equal to the double compare. 
    bool diameter(double compare){
        //get widths of intervals and store them in mpfr variable
        mpfi_diam_abs(maxDiameterX,  x);
        mpfi_diam_abs(maxDiameterY,  y);
        mpfi_diam_abs(maxDiameterZ, z);
        
        //get max values between each of the mfpr widths
        mpfr_max((mpfr_ptr)intermediateMax1, maxDiameterX, maxDiameterY, MPFR_RNDN);
        mpfr_max((mpfr_ptr) intermediateMax2, intermediateMax1, maxDiameterZ, MPFR_RNDN);

        //compare mpfr max with the compare double value and return. if intermediateMax2 compare == 1, if they are equal, 
        //compare == 0, else compare == -1 so we need to correct for this in an if-else statement. 
        int truth = mpfr_cmp_d(intermediateMax2, compare);
        if(truth <= 0){
            truth = true;
        }
        else{
            truth = false;
        }

        return truth;
    }

    //clear intervals after program is used. (All mpfi_t variables should be cleared according to docs to free memory)
    void Clear() {
        mpfi_clear((mpfi_ptr)x);
        mpfi_clear((mpfi_ptr)y);
        mpfi_clear((mpfi_ptr)z);

        mpfr_clear(maxDiameterX);
        mpfr_clear(maxDiameterY);
        mpfr_clear(maxDiameterZ);

        mpfr_clear(intermediateMax1);
        mpfr_clear(intermediateMax2);

    }

};






/*
function that defines poly 1. It performs interval arithmetic on the given x, y, z intervals
then stores the calculated interval in ROP (resultant output). Each interval operation in the 
library is defined as a function which stores the output of the operation in a new, separately provided mpfi_t interval
(designated ROP) to the function. Thus, each time I perform an operation in the following function, I first initialize
a new, intermediate interval that stores the intermediate result. 

This is an example of a sphere
*/

void f(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){
    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);


    mpfi_t zSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) zSquared);
    //square z
    mpfi_sqr(zSquared, z);


    //Adding x^2 and y^2
    mpfi_t firstAddition;
    mpfi_init((mpfi_ptr) firstAddition);
    mpfi_add ((mpfi_ptr) firstAddition, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);

    //Adding firstAddition and z^2
    mpfi_t secondAddition;
    mpfi_init((mpfi_ptr) secondAddition);
    mpfi_add ((mpfi_ptr) secondAddition, (mpfi_ptr) firstAddition, (mpfi_ptr) zSquared);

    //third and final output. Assign the output to ROP 
    


    mpfi_sub_d(ROP, secondAddition, 1);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) zSquared);
    mpfi_clear((mpfi_ptr) firstAddition);
    mpfi_clear((mpfi_ptr) secondAddition);
}

/*
function that defines poly 2. It performs interval arithmetic on the given x, y, z intervals
then stores the calculated interval in ROP (resultant output). Each interval operation in the 
library is defined as a function which stores the output of the operation in a new, separately provided mpfi_t interval
(designated ROP) to the function. Thus, each time I perform an operation in the following function, I first initialize
a new, intermediate interval that stores the intermediate result. 

This is an example of an elliptic paraboloid. 
*/

void g(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){
    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);

    //Adding x^2 and y^2
    mpfi_t firstAddition;
    mpfi_init((mpfi_ptr) firstAddition);
    mpfi_add ((mpfi_ptr) firstAddition, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);

    //Subtracting z from firstAddition 
    mpfi_t secondSubtraction;
    mpfi_init((mpfi_ptr) secondSubtraction);
    mpfi_sub ((mpfi_ptr) secondSubtraction, (mpfi_ptr) firstAddition, (mpfi_ptr) z);

    //third and final output. Assign the output to ROP 
    
    mpfi_sub_d(ROP, secondSubtraction, 1);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) firstAddition);
    mpfi_clear((mpfi_ptr) secondSubtraction);
}






// hyperboloidOfOneSheet function in the form of x^2/a^2 + y^2/b^2 - z^2/c^2 - 1. Makesure to clear ROP, x, y, z somewhere else. Stores result of evaluated function in ROP
void hyperboloidOfOneSheet(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){

    //initializing constants. 
    double a = 1;
    double b = 0.75; 
    double c = 1;

    a = pow(a, 2);
    b = pow(b, 2);
    c = pow(c, 2);


    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);


    mpfi_t zSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) zSquared);
    //square z
    mpfi_sqr(zSquared, z);


    //Dividing each of the squared interavls by their respective constants (check to make sure can pass the same argument)
    mpfi_div_d(xSquared, xSquared, a);
    mpfi_div_d(ySquared, ySquared, b);
    mpfi_div_d(zSquared, zSquared, c);


    //Adding x^2 and y^2
    mpfi_t firstAddition;
    mpfi_init((mpfi_ptr) firstAddition);
    mpfi_add ((mpfi_ptr) firstAddition, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);

    //subtracting z^2/c from firstAddition 
    mpfi_t secondSubtraction;
    mpfi_init((mpfi_ptr) secondSubtraction);
    mpfi_sub ((mpfi_ptr) secondSubtraction, (mpfi_ptr) firstAddition, (mpfi_ptr) zSquared);

    //third and final output. Assign the output to ROP 
    


    mpfi_sub_d(ROP, secondSubtraction, 1);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) zSquared);
    mpfi_clear((mpfi_ptr) firstAddition);
    mpfi_clear((mpfi_ptr) secondSubtraction);
    
}

// hyperboloidOfTwoSheets function in the form of x^2/a^2 + y^2/b^2 - z^2/c^2 + 1. Makesure to clear ROP, x, y, z somewhere else. Stores result of evaluated function in ROP
void hyperboloidOfTwoSheets(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){

    //initializing constants. 
    double a = 1;
    double b = 2; 
    double c = 0.5;

    a = pow(a, 2);
    b = pow(b, 2);
    c = pow(c, 2);


    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);


    mpfi_t zSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) zSquared);
    //square z
    mpfi_sqr(zSquared, z);


    //Dividing each of the squared interavls by their respective constants (check to make sure can pass the same argument)
    mpfi_div_d(xSquared, xSquared, a);
    mpfi_div_d(ySquared, ySquared, b);
    mpfi_div_d(zSquared, zSquared, c);


    //Adding x^2 and y^2
    mpfi_t firstAddition;
    mpfi_init((mpfi_ptr) firstAddition);
    mpfi_add ((mpfi_ptr) firstAddition, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);

    //subtracting z^2/c from firstAddition 
    mpfi_t secondSubtraction;
    mpfi_init((mpfi_ptr) secondSubtraction);
    mpfi_sub ((mpfi_ptr) secondSubtraction, (mpfi_ptr) firstAddition, (mpfi_ptr) zSquared);

    //third and final output. Assign the output to ROP 
    
    mpfi_add_d(ROP, secondSubtraction, 1);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) zSquared);
    mpfi_clear((mpfi_ptr) firstAddition);
    mpfi_clear((mpfi_ptr) secondSubtraction);
    
}


//calculates equation for hyperplane aX + bY + cZ + d = 0. Make sure to clear ROP, x, y, z somewhere else
void hyperplane(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){
    double a = 0; 
    double b = 0 ;
    double c = -1;
    double d = 0.5;


    mpfi_t xScaled;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xScaled);
    //square x
    mpfi_mul_d(xScaled, x, a);

    mpfi_t yScaled;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) yScaled);
    //square y
    mpfi_mul_d(yScaled, y, b);


    mpfi_t zScaled;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) zScaled);
    //square z
    mpfi_mul_d(zScaled, z, c);


    //Adding x^2 and y^2
    mpfi_t firstAddition;
    mpfi_init((mpfi_ptr) firstAddition);
    mpfi_add ((mpfi_ptr) firstAddition, (mpfi_ptr) xScaled, (mpfi_ptr) yScaled);

    //Adding firstAddition and z^2
    mpfi_t secondAddition;
    mpfi_init((mpfi_ptr) secondAddition);
    mpfi_add ((mpfi_ptr) secondAddition, (mpfi_ptr) firstAddition, (mpfi_ptr) zScaled);

    //third and final output. Assign the output to ROP 
    


    mpfi_add_d(ROP, secondAddition, d);

    //clearing variables
    mpfi_clear((mpfi_ptr) xScaled);
    mpfi_clear((mpfi_ptr) yScaled);
    mpfi_clear((mpfi_ptr) zScaled);
    mpfi_clear((mpfi_ptr) firstAddition);
    mpfi_clear((mpfi_ptr) secondAddition);

}



// elipsoid function in the form of x^2/a^2 + y^2/b^2 + z^2/c^2 - 1 = 0. Makesure to clear ROP, x, y, z somewhere else. Stores result of evaluated function in ROP
void ellipsoid(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){

    //initializing constants. 
    double a = 1;
    double b = 2; 
    double c = 1;

    a = pow(a, 2);
    b = pow(b, 2);
    c = pow(c, 2);

    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);


    mpfi_t zSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) zSquared);
    //square z
    mpfi_sqr(zSquared, z);


    //Dividing each of the squared interavls by their respective constants (check to make sure can pass the same argument)
    mpfi_div_d(xSquared, xSquared, a);
    mpfi_div_d(ySquared, ySquared, b);
    mpfi_div_d(zSquared, zSquared, c);


    //Adding x^2 and y^2
    mpfi_t firstAddition;
    mpfi_init((mpfi_ptr) firstAddition);
    mpfi_add ((mpfi_ptr) firstAddition, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);

    //subtracting z^2/c from firstAddition 
    mpfi_t secondAddition;
    mpfi_init((mpfi_ptr) secondAddition);
    mpfi_add ((mpfi_ptr) secondAddition, (mpfi_ptr) firstAddition, (mpfi_ptr) zSquared);

    //third and final output. Assign the output to ROP 
    


    mpfi_sub_d(ROP, secondAddition, 1);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) zSquared);
    mpfi_clear((mpfi_ptr) firstAddition);
    mpfi_clear((mpfi_ptr) secondAddition);
    
}

// elliptic Cone function in the form of x^2/a^2 + y^2/b^2 - z^2/c^2= 0. Makesure to clear ROP, x, y, z somewhere else. Stores result of evaluated function in ROP
void ellipticCone(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){

    //initializing constants. 
    double a = 2;
    double b = 3; 
    double c = 2;

    a = pow(a, 2);
    b = pow(b, 2);
    c = pow(c, 2);

    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);


    mpfi_t zSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) zSquared);
    //square z
    mpfi_sqr(zSquared, z);


    //Dividing each of the squared interavls by their respective constants (check to make sure can pass the same argument)
    mpfi_div_d(xSquared, xSquared, a);
    mpfi_div_d(ySquared, ySquared, b);
    mpfi_div_d(zSquared, zSquared, c);


    //Adding x^2 and y^2
    mpfi_t firstAddition;
    mpfi_init((mpfi_ptr) firstAddition);
    mpfi_add ((mpfi_ptr) firstAddition, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);

    //subtracting z^2/c^2 from firstAddition. Storing output in ROP
    
    mpfi_sub((mpfi_ptr) ROP, (mpfi_ptr) firstAddition, (mpfi_ptr) zSquared);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) zSquared);
    mpfi_clear((mpfi_ptr) firstAddition);
    
}

/*
function that defines hyperbolic paraboloid. It performs interval arithmetic on the given x, y, z intervals
then stores the calculated interval in ROP (resultant output). The equation for a hyperbolic paraboloid is given by x^2/a^2 - y^2/b^2 = z/c
*/

void hyperbolicParaboloid(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){
    //constants a,b,c
    double a = 1; 
    double b = 1; 
    double c = 1; 

    a = pow(a,2);
    b = pow(b,2);



    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_div_d(xSquared, xSquared, a);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);

    mpfi_div_d(ySquared, ySquared, b);

    //dividing z by c

    mpfi_t zDividedByC;
    mpfi_init(zDividedByC);

    mpfi_div_d(zDividedByC, z, c);


    //subtracting y^2 from x^2
    mpfi_t firstSubtraction;
    mpfi_init((mpfi_ptr) firstSubtraction);
    mpfi_sub ((mpfi_ptr) firstSubtraction, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);

    //Subtracting z from firstSubtraction and then storing result in ROP 
    
    mpfi_sub ((mpfi_ptr) ROP, (mpfi_ptr) firstSubtraction, (mpfi_ptr) zDividedByC);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear(zDividedByC);
    mpfi_clear((mpfi_ptr) firstSubtraction);
}

/*
function that defines an elliptic paraboloid of the form x^2/a^2 + y^2/b^2 - z = 0. It performs interval arithmetic on the given x, y, z intervals
then stores the calculated interval in ROP (resultant output). 
*/

void ellipticParabolid(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){
    //constants a,b,c
    double a = 2; 
    double b = 2; 
    

    a = pow(a,2);
    b = pow(b,2);



    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_div_d(xSquared, xSquared, a);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);

    mpfi_div_d(ySquared, ySquared, b);

    //Adding x^2 and y^2
    mpfi_t firstSubtraction;
    mpfi_init((mpfi_ptr) firstSubtraction);
    mpfi_add ((mpfi_ptr) firstSubtraction, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);

    //Subtracting z from firstSubtraction and then storing result in ROP 
    
    mpfi_sub ((mpfi_ptr) ROP, (mpfi_ptr) firstSubtraction, (mpfi_ptr) z);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) firstSubtraction);
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
/**
 * The below functions of conic sections return boolean values 
*/
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////



/***
 * This is a function calculating whether for a given x,y, and z interval, the equation for an elliptic cylinder is defined. It does so by returning a boolean value
 * corresponding to 1 if the cylinder is defined, and 0 if not. Note that unlike the other functions defined in this program, this is not void but
 * returns a boolean value to account for the fact that it must test if the z interval is within the given z axis bounds of the cylinder equation.   
 * z axis bounds are optional. General equation of a cylinder is x^2 + y^2 = r
 * 
 * 
*/

bool ellipticCylinder(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){

    //scaling x^2 and y^2 
    double a = 1;
    double b = 0.8;

    a = pow(a, 2);
    b = pow(b, 2);

    //initializing constants. 
    //lower and upper bounds in z axis
    double z1 = 0;
    double z2 = 1; 

    //radius
    double r = 0.5;


    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_div_d(xSquared, xSquared, a);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);

    mpfi_div_d(ySquared, ySquared, b);


        
    //Adding x^2 and y^2
    mpfi_t firstAddition;
    mpfi_init((mpfi_ptr) firstAddition);
    mpfi_add ((mpfi_ptr) firstAddition, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);


    //final output. Assign the output to ROP. Subtracting radius from first addition
    
    mpfi_sub_d(ROP, firstAddition, r);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) firstAddition);

    //check if ROP has 0 contained. If so, then check bounds on z1 and z2
    if(mpfi_has_zero(ROP)){

        //Getting the left and right endpoints of the z interval
        mpfr_t leftEndpoint; 
        mpfr_t rightEndpoint; 
        mpfr_init(leftEndpoint);
        mpfr_init(rightEndpoint);
        mpfi_get_left(leftEndpoint, z);
        mpfi_get_right(rightEndpoint, z);

        //comparing the endpoints to z1 and z2. -1 if endpoint is smaller, 0 if equal, 1 if endpoint is larger. 
        int compareZ2 = mpfr_cmp_d(leftEndpoint, z2);
        int compareZ1 = mpfr_cmp_d(rightEndpoint, z1);
        


        //case 1: if left endpoint <= z2 (compareZ2 <= 0) and right endpoint >= z1 (compareZ1 >= 0)
        if(( compareZ2 <= 0) && (compareZ1 >= 0) ){
            mpfr_clear(leftEndpoint);
            mpfr_clear(rightEndpoint);
            return true;
        }
        else{ //case 2: the z interval does not contain any points between endpoints z1 and z2
            mpfr_clear(leftEndpoint);
            mpfr_clear(rightEndpoint);
            return false;
        }
    }
    else{
        return false;
    }

}

/***
 * This is a function calculating whether for a given x,y, and z interval, the equation for an parabolic cylinder is defined. It does so by returning a boolean value
 * corresponding to 1 if the parabolic cylinder is defined, and 0 if not. Note that unlike the other functions defined in this program, this is not void but
 * returns a boolean value to account for the fact that it must test if the z interval is within the given z axis bounds of the cylinder equation.   
 * z axis bounds are optional. General equation of a parabolic cylinder is y = ax^2
*/

bool parabolicCylinder(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){

    //scaling x^2 and y^2 
    double a = 1;

    //initializing constants. 
    //lower and upper bounds in z axis
    double z1 = -1;
    double z2 = 1; 

    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_mul_d(xSquared, xSquared, a);
        
    //subtracting ax^2 from y^2 and storing it in ROP
    
    mpfi_sub ((mpfi_ptr) ROP, (mpfi_ptr) y, xSquared);


    

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);

    //check if ROP has 0 contained. If so, then check bounds on z1 and z2
    if(mpfi_has_zero(ROP)){

        //Getting the left and right endpoints of the z interval
        mpfr_t leftEndpoint; 
        mpfr_t rightEndpoint; 
        mpfr_init(leftEndpoint);
        mpfr_init(rightEndpoint);
        mpfi_get_left(leftEndpoint, z);
        mpfi_get_right(rightEndpoint, z);

        //comparing the endpoints to z1 and z2. -1 if endpoint is smaller, 0 if equal, 1 if endpoint is larger. 
        int compareZ2 = mpfr_cmp_d(leftEndpoint, z2);
        int compareZ1 = mpfr_cmp_d(rightEndpoint, z1);
        


        //case 1: if left endpoint <= z2 (compareZ2 <= 0) and right endpoint >= z1 (compareZ1 >= 0)
        if(( compareZ2 <= 0) && (compareZ1 >= 0) ){
            mpfr_clear(leftEndpoint);
            mpfr_clear(rightEndpoint);
            return true;
        }
        else{ //case 2: the z interval does not contain any points between endpoints z1 and z2
            mpfr_clear(leftEndpoint);
            mpfr_clear(rightEndpoint);
            return false;
        }
    }
    else{
        return false;
    }

}


/***
 * This is a function calculating whether for a given x,y, and z interval, the equation for a hyperbolic cylinder is defined. It does so by returning a boolean value
 * corresponding to 1 if the cylinder is defined, and 0 if not. Note that unlike the other functions defined in this program, this is not void but
 * returns a boolean value to account for the fact that it must test if the z interval is within the given z axis bounds of the cylinder equation.   
 * z axis bounds are optional. General equation of a hyperbolic cylinder is x^2/a^2 - y^2/b^2 = 1
*/

bool hyperbolicCylinder(mpfi_t ROP,mpfi_t x, mpfi_t y, mpfi_t z){

    //scaling x^2 and y^2 
    double a = 0.5;
    double b = 0.5;

    a = pow(a, 2);
    b = pow(b, 2);

    //initializing constants. 
    //lower and upper bounds in z axis
    double z1 = -1;
    double z2 = 1; 

    mpfi_t xSquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) xSquared);
    //square x
    mpfi_sqr(xSquared, x);

    mpfi_div_d(xSquared, xSquared, a);

    mpfi_t ySquared;
    // Initialize the interval variable
    mpfi_init((mpfi_ptr) ySquared);
    //square y
    mpfi_sqr(ySquared, y);

    mpfi_div_d(ySquared, ySquared, b);


        
    //subtracting y^2 from x^2
    mpfi_t firstSubtraction;
    mpfi_init((mpfi_ptr) firstSubtraction);
    mpfi_sub ((mpfi_ptr) firstSubtraction, (mpfi_ptr) xSquared, (mpfi_ptr) ySquared);


    //final output. Assign the output to ROP. Subtracting radius from first addition
    
    mpfi_sub_d(ROP, firstSubtraction, 1);

    //clearing variables
    mpfi_clear((mpfi_ptr) xSquared);
    mpfi_clear((mpfi_ptr) ySquared);
    mpfi_clear((mpfi_ptr) firstSubtraction);

    //check if ROP has 0 contained. If so, then check bounds on z1 and z2
    if(mpfi_has_zero(ROP)){

        //Getting the left and right endpoints of the z interval
        mpfr_t leftEndpoint; 
        mpfr_t rightEndpoint; 
        mpfr_init(leftEndpoint);
        mpfr_init(rightEndpoint);
        mpfi_get_left(leftEndpoint, z);
        mpfi_get_right(rightEndpoint, z);

        //comparing the endpoints to z1 and z2. -1 if endpoint is smaller, 0 if equal, 1 if endpoint is larger. 
        int compareZ2 = mpfr_cmp_d(leftEndpoint, z2);
        int compareZ1 = mpfr_cmp_d(rightEndpoint, z1);
        


        //case 1: if left endpoint <= z2 (compareZ2 <= 0) and right endpoint >= z1 (compareZ1 >= 0)
        if(( compareZ2 <= 0) && (compareZ1 >= 0) ){
            mpfr_clear(leftEndpoint);
            mpfr_clear(rightEndpoint);
            return true;
        }
        else{ //case 2: the z interval does not contain any points between endpoints z1 and z2
            mpfr_clear(leftEndpoint);
            mpfr_clear(rightEndpoint);
            return false;
        }
    }
    else{
        return false;
    }

}


/*********
 * 
 * 
 *
 * This code chunk below is for the manual expression functionality. It accepts an arbitrary sum/difference of monomials
 * 
 * 
 */

//This seperates an expression into operators and monomials

void separateExpression(std::string expression, std::vector<std::string> &variables, std::vector<std::string> &arr) {
    std::stringstream ss(expression);
    std::string token;
    bool expectOperator = false;

    // the >> operator is usually right shift, but in stream contexts can act as "read until next delimiter (in this program, it is 
    //white space). Then stores this in token. "
    while (ss >> token) {

        if (token == "+" || token == "-") {
            //next token should 
            if (!expectOperator) {
            std::cerr << "Wrong input: Unexpected operator at this position (should've been an algebraic expression).\n";
            return;
            }
            //If the token is an operator, add it to the arr vector


            arr.push_back(token);

            expectOperator = false;
        } else {

            if (expectOperator) {
            std::cerr << "Wrong input: Unexpected expression at this position (should've been an operator).\n";
            return;
            }

            // If the token is a monomial, add it to the variables vector
            variables.push_back(token);
            expectOperator = true;
        }
    }
}

//extracts the power a particular variable is raised to. 
int extractPower(std::string str, char variable) {

    //use size_t because that is how sizes are defined in c++ and also it's default return type of these functions anyway
    //find the position of the variable
    size_t pos = str.find(variable);
    if (pos == std::string::npos) {
        return -1;
    }  
    
    //find the position of '^' after the variable
    size_t pow_pos = str.find('^', pos);
    if (pow_pos == std::string::npos) {
         return -1;
        
    
    }  
    //find the position of the next variable or the end of the string
    size_t next_var_pos = str.find_first_of("xyz", pow_pos + 1);
    
    //extract the power substring. substr(starting index, length) 
    std::string power_str = str.substr(pow_pos + 1, next_var_pos - pow_pos - 1);
    
    //convert the substring to an integer
    return std::stoi(power_str);
}

//raises an interval to an arbitrary power

void raiseToPower(mpfi_t result, mpfi_t inputInterval, int exponentNumerator, int exponentDenominator) {
    mpfi_t temp;
    mpfi_init(temp);
    //let p = exponent
    mpfi_set(temp, inputInterval);


    //loop through power - 1 times since we've already set temp equal to the first power
    for(int i = 0; i < exponentNumerator -1; i++){
        mpfi_mul(temp, temp, inputInterval);
    }


    //checking if the interval is raised to an even power. 
    //if raised to an even power, then check if contains negative values. If 0 is contained, then assumed to have both negative and 
    //positive values, if no 0, than just return interval
    if(exponentNumerator %2 == 0){
        //if has 0, set left endpoint to 0, and keep right endpoint. If not, do nothing and leave interval since there should be 
        //no negative values in that case. 
        if(mpfi_has_zero(temp)){
            mpfr_t rightEndpoint;
            mpfr_init(rightEndpoint);
            mpfi_get_right(rightEndpoint, temp);
            //rounding away from 0 (encompoass larger margin of error always)
            double rightEndpointInDoubleFormat = mpfr_get_d(rightEndpoint, MPFR_RNDU);
            mpfi_interv_d(temp, 0.0, rightEndpointInDoubleFormat);
            mpfr_clear(rightEndpoint);
        }
        
    }

   

    mpfi_set(result, temp);
    mpfi_clear(temp);
}


//returns what case the expresion is. We use this case to determine how to evaluate an expression on a given interval

// the following case numbers correspond to each possibility: 

// x => 1, y => 2, x, xy => 3 , z => 4, xz => 5, yz => 6, xyz => 7, no variables = 0

int caseNumber(std::string expression){

    int code = 0;

    // Check for the presence of 'x', 'y', 'z'
    if (expression.find('x') != std::string::npos) code |= 1; // Set bit 0
    if (expression.find('y') != std::string::npos) code |= 2; // Set bit 1
    if (expression.find('z') != std::string::npos) code |= 4; // Set bit 2

    return code;



}


//without parenthesis 

// sets resultOfExpression (assumed to be initialized but not set interval) to the result of the expression when evaluated on the 
//given x,y,z intervals. 

//can't handle fractional powers. 

void eval(mpfi_t resultOfExpression, std::string expression, mpfi_t xInterval, mpfi_t yInterval, mpfi_t zInterval){

    //check if the there is a constant at the start of the expression
    double constant;
    if(isdigit(expression[0])){
        //extract the constant number
        size_t next_var_pos = expression.find_first_of("xyz",  0);
    
        //extract the power substring. substr(starting index, length) 
        constant = std::stod((expression.substr(0, next_var_pos))) ;
        
    }

    //std::cout << "This expression :" << expression << "\n\n\n";


    int code = caseNumber(expression);

    
    switch (code ) {

        //case x^something
        case 1: {
            raiseToPower(resultOfExpression, xInterval, extractPower(expression,'x'), 1);

            if(isdigit(expression[0])){
                mpfi_mul_d(resultOfExpression,resultOfExpression, constant);   
            }

            return;
        }

        // Case 'y^something' 
        case 2: {
            raiseToPower(resultOfExpression, yInterval, extractPower(expression,'y'), 1);

            if(isdigit(expression[0])){
                mpfi_mul_d(resultOfExpression,resultOfExpression, constant);   
            }
            return ; 
        }
        // Case 'z^something'
        case 4: {

            raiseToPower(resultOfExpression, zInterval, extractPower(expression,'z'), 1);

            if(isdigit(expression[0])){
                mpfi_mul_d(resultOfExpression,resultOfExpression, constant);   
            }

            
            return ; 
        }
        // Case 'x^something y^something'
        case 3: {


            mpfi_t xRaisedToPower;
            mpfi_init(xRaisedToPower);
            mpfi_t yRaisedToPower;
            mpfi_init(yRaisedToPower);

            raiseToPower(xRaisedToPower, xInterval, extractPower(expression, 'x'), 1);
            raiseToPower(yRaisedToPower, yInterval, extractPower(expression, 'y'), 1);

            mpfi_mul(resultOfExpression, xRaisedToPower, yRaisedToPower);
            if(isdigit(expression[0])){
                mpfi_mul_d(resultOfExpression,resultOfExpression, constant);   
            }

            mpfi_clear(xRaisedToPower);
            mpfi_clear(yRaisedToPower);
            return ; 
        }
         // Case 'x^something z^something'
        case 5: {


            mpfi_t xRaisedToPower;
            mpfi_init(xRaisedToPower);
            mpfi_t zRaisedToPower;
            mpfi_init(zRaisedToPower);


            raiseToPower(xRaisedToPower, xInterval, extractPower(expression, 'x'), 1);
            raiseToPower(zRaisedToPower, zInterval, extractPower(expression, 'z'), 1);

            mpfi_mul(resultOfExpression, xRaisedToPower, zRaisedToPower);
            if(isdigit(expression[0])){
                mpfi_mul_d(resultOfExpression,resultOfExpression, constant);   
            }

            mpfi_clear(xRaisedToPower);
            mpfi_clear(zRaisedToPower);
            return;
        }

        // Case 'y^something z^something'
        case 6: {

            mpfi_t yRaisedToPower;
            mpfi_init(yRaisedToPower);
            mpfi_t zRaisedToPower;
            mpfi_init(zRaisedToPower);


            raiseToPower(yRaisedToPower, yInterval, extractPower(expression, 'y'), 1);
            raiseToPower(zRaisedToPower, zInterval, extractPower(expression, 'z'), 1);

            mpfi_mul(resultOfExpression, yRaisedToPower, zRaisedToPower);
            if(isdigit(expression[0])){
                mpfi_mul_d(resultOfExpression,resultOfExpression, constant);   
            }

            mpfi_clear(yRaisedToPower);
            mpfi_clear(zRaisedToPower);
            
            return ; 
        }

        // Case 'x^something y^something z^something'
        case 7: {


            mpfi_t xRaisedToPower;
            mpfi_init(xRaisedToPower);
            mpfi_t yRaisedToPower;
            mpfi_init(yRaisedToPower);
            mpfi_t zRaisedToPower;
            mpfi_init(zRaisedToPower);


            raiseToPower(xRaisedToPower, xInterval, extractPower(expression, 'x'), 1);
            raiseToPower(yRaisedToPower, yInterval, extractPower(expression, 'y'), 1);
            raiseToPower(zRaisedToPower, zInterval, extractPower(expression, 'z'), 1);


            mpfi_mul(resultOfExpression, xRaisedToPower, yRaisedToPower);
            mpfi_mul(resultOfExpression, resultOfExpression, zRaisedToPower);
            if(isdigit(expression[0])){
                mpfi_mul_d(resultOfExpression,resultOfExpression, constant);   
            }

            mpfi_clear(xRaisedToPower);
            mpfi_clear(zRaisedToPower);
            mpfi_clear(yRaisedToPower);
            return;
        }

        //case where just a constant
        case 0:
        {
            //set interval to the value of a double
            mpfi_interv_d(resultOfExpression,std::stod(expression), std::stod(expression));   
            
        
            return ; 
        }
    }
}

void abstractExpressionCalculator(std::string expression, mpfi_t x, mpfi_t y, mpfi_t z, mpfi_t finalAnswer){


    

    //first parse expression and fill the operations vector (arr) and the variables vector
    std::vector<std::string>  arr ; 

    std::vector<std::string> variables;

    separateExpression(expression, variables, arr);

    // for(const auto &i: arr){
    //     std::cout << i << " ";
    // }

    // std::cout << "\n\n";

    // for(const auto &i: variables){
    //     std::cout << i << " ";
    // }

    std::cout << "\n\n";

    //start iterating through the expression 
    mpfi_t tempAddOrSub;

    mpfi_init(tempAddOrSub);

    //evaluate first part of expression
    eval(finalAnswer, variables[0], x, y, z);

    mpfi_out_str(stdout, 10, 10, finalAnswer);
    std::cout << "\n";


    //use eval function for each variable
    for(int i = 1; i < variables.size(); i++){

        


        eval(tempAddOrSub, variables[i], x, y, z);



        //adding and subtracting: 
        if(arr[i-1] == "+"){
            mpfi_add(finalAnswer, finalAnswer, tempAddOrSub);
        }
        else{
            mpfi_sub(finalAnswer, finalAnswer, tempAddOrSub);

        }       
    }

    
    

    mpfi_clear(tempAddOrSub);

}


/*********
 * 
 * 
 *
 * This code chunk above is for the manual expression functionality
 * 
 * 
 */


//function that checks if 0 is defined for both f and g
bool hasZero(mpfi_t x, mpfi_t y, mpfi_t z) {

    // initializing intervals to store f and g, and then computing f and g with the given x,y,z intervals. 
    mpfi_t resultOfF;
    mpfi_init((mpfi_ptr) resultOfF);

    mpfi_t resultOfg;
    mpfi_init((mpfi_ptr) resultOfg);


    //this is an optional expression in case a manual expression would like to be used. 

    //dogan example 1
    std::string expression = "x^4 + 2x^2y^2 + y^4 - 2x^2 - 2y^2 + 1 - z^1";

    //dogan example 2
    std::string doganExperiment4a = "x^2 + y^2 - z^2 - 2";
    std::string doganExperiment4b = "x^2 - y^2 + z^2 - 1";


    //This is the part of the code that you change to calculate for various function intersections. 

    //f(resultOfF, x, y, z);

    

    //g(resultOfg, x, y, z);

    //For visualization

    //f(resultOfF, x, y, z);
    abstractExpressionCalculator(doganExperiment4a, x, y, z, resultOfg);
    abstractExpressionCalculator(doganExperiment4b, x, y, z, resultOfF);

    // hyperplane(resultOfF, x, y, z);
    

    //testing to make sure each calculation is the same: 

    // mpfr_t rightEndpoint;
    // mpfr_init(rightEndpoint);
    // mpfr_t rightEndpoint2;
    // mpfr_init(rightEndpoint2);
    // mpfr_t leftEndpoint;
    // mpfr_init(leftEndpoint);
    // mpfr_t leftEndpoint2;
    // mpfr_init(leftEndpoint2);



    // mpfi_get_right(rightEndpoint, resultOfF);
    // mpfi_get_right(rightEndpoint2, resultOfg);
    // mpfi_get_left(leftEndpoint, resultOfF);
    // mpfi_get_left(leftEndpoint2, resultOfg);

    // //rounding away from 0 (encompoass larger margin of error always)
    // double rightEndpointInDoubleFormat = mpfr_get_d(rightEndpoint, MPFR_RNDU);
    // double right2EndpointInDoubleFormat = mpfr_get_d(rightEndpoint2, MPFR_RNDU);
    // double leftEndpointInDoubleFormat = mpfr_get_d(leftEndpoint, MPFR_RNDU);
    // double left2EndpointInDoubleFormat = mpfr_get_d(leftEndpoint2, MPFR_RNDU);

    // if(rightEndpointInDoubleFormat != right2EndpointInDoubleFormat || leftEndpointInDoubleFormat != left2EndpointInDoubleFormat){
    //     std::cout << "PROBLEM: Here are the input intervals, then after are the calculations\n";
    //     std::cout << " These next three are the x,y,z intervals\n";

    // mpfi_out_str(stdout, 10, 10, x);
    // std::cout << "\n";

    //  mpfi_out_str(stdout, 10, 10, y);
    // std::cout << "\n";


    // mpfi_out_str(stdout, 10, 10, z);
    // std::cout << "\n";

    // mpfi_out_str(stdout, 10, 10, resultOfF);
    // std::cout << "\n";


    // mpfi_out_str(stdout, 10, 10, resultOfg);
    // std::cout << "\n";
    // }


    // mpfr_clear(rightEndpoint);
    // mpfr_clear(rightEndpoint2);
    // mpfr_clear(leftEndpoint);
    // mpfr_clear(leftEndpoint2);


    //end test
    


    



    //checking if f and g both have zero in them. Returning true if both have zero. 
    if(mpfi_has_zero(resultOfF) && mpfi_has_zero(resultOfg)){
        mpfi_clear((mpfi_ptr) resultOfF);
        mpfi_clear((mpfi_ptr) resultOfg);
        return true;
    }
    else{
        mpfi_clear((mpfi_ptr) resultOfF);
        mpfi_clear((mpfi_ptr) resultOfg);
        return false;
    }
}



int main(){




    

    //define interval box that I will be computing within 
    intervalClass initialInterval = intervalClass(-3, 3, -3, 3, -3, 3);


    //define error
    double maxError = 0.5 / (pow(2.0, 4));
    




    //initialize queue of all intervals 
    std::queue<intervalClass> intervalQueue;
    //update with the first box
    intervalQueue.push(initialInterval);

    //initialize empty answer vector of intervals to contain where f and g are both equal to 0
    std::vector<intervalClass> answerIntervals;

    //while queue is nonempty, run code. To pop from queue, the width of its intervals must be less than or equal to the maxError

    while(!intervalQueue.empty()){
        intervalClass currBox = intervalQueue.front();
        intervalQueue.pop();

        

    
    //inside while
    //pop front of queue, set equal to a variable 
    //check if f and g are equal to 0 in this box
    if(hasZero(currBox.x, currBox.y, currBox.z)){
    //check if current box has all dimensions less than desired error:
    if(currBox.diameter(maxError)){
        //if yes, then append to answer array
        answerIntervals.push_back(currBox);

    }
    else{
        //else append newly created boxes (8 total) to queue
        //first initialize arrays
        mpfi_t newX[2];
        mpfi_t newY[2];
        mpfi_t newZ[2];

        //initializing mpfi_t arrays before placing values into them
        for(int i = 0; i<2; i++){
            mpfi_init(newX[i]);
            mpfi_init(newY[i]);
            mpfi_init(newZ[i]);
        }



        mpfi_t x; 
        mpfi_init(x);
        mpfi_set(x, currBox.x);

        mpfi_t y; 
        mpfi_init(y);
        mpfi_set(y, currBox.y);

        mpfi_t z; 
        mpfi_init(z);
        mpfi_set(z, currBox.z);

        //splitting box into 8 now boxes, so splitting x,y,z into 2 each so we can get 2 * 2 * 2 = 8 different combinations of boxes
        mpfi_bisect((mpfi_ptr)newX[0], (mpfi_ptr) newX[1], x);
        mpfi_bisect((mpfi_ptr)newY[0], (mpfi_ptr) newY[1], y);
        mpfi_bisect((mpfi_ptr)newZ[0], (mpfi_ptr) newZ[1], z);

        //triple for loop to queue all new elements. 
        for(int i = 0; i <2 ; i++){
            for(int j = 0; j < 2; j++){
                for(int k = 0; k <2; k++){
                    intervalQueue.push(intervalClass(newX[i], newY[j], newZ[k]));
                }
            }
        }

        //clearing out variables that aren't needed
        currBox.Clear();
        for(int i = 0; i<2; i++){
            mpfi_clear((mpfi_ptr) newX[i]);
            mpfi_clear((mpfi_ptr) newY[i]);
            mpfi_clear((mpfi_ptr)newZ[i]);
        }
       
    }
    
    }
    else{

        //if no, continue to next iteration of while loop and clear the box
        currBox.Clear();
        continue;
    }
    }
    //outside while
    //print/postprocess ans array, make sure to clear boxes in ans array too. 


    //Write to file

    FILE *fptr;

    // Open a file in append mode
    fptr = fopen("/Users/sohumbelagode/Desktop/Yap- Research/initialProject.txt", "w");

    // Append some text to the file
    fprintf(fptr, "The boxes: \n");
    // Close the file
    fclose(fptr);

    fptr = fopen("/Users/sohumbelagode/Desktop/Yap- Research/initialProject.txt", "a");


    long int counter = 0;
    for (std::vector<intervalClass>::iterator it = answerIntervals.begin(); it != answerIntervals.end(); ++it) {
        intervalClass box = *it; 
        counter++;
        
        fprintf(fptr, "Box ");
        fprintf(fptr,"%ld", counter); 
        fprintf(fptr, " : \n");

        std::cout << "Box " << counter + 1 << ": "<< std::endl;
        box.p();
        box.p(fptr);

        std::cout << std::endl;
        fprintf(fptr, "\n");
        box.Clear();
        
    }


    std::cout << "Total Number of Boxes: " << counter << " Error: "<< maxError<<std::endl;

    fprintf(fptr, "Total Number of Boxes: ");
    fprintf(fptr,"%ld", counter); 
    fprintf(fptr, " Error: ");
    fprintf(fptr,"%lf", maxError); 





    

    fclose(fptr);
    return 0;
} 


