
#include <iostream>
//import library for computing interval math
#include "mpfi.h"
#include "mpfi_io.h"
#include <queue>

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
    }


    

    //print method

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
    }

    //gets max diameter and then checks if it is less than or equal to the double compare. 
    bool diameter(double compare){
        //get widths of intervals and store them in mpfr variable
        mpfi_diam_abs(maxDiameterX, (mpfi_ptr) x);
        mpfi_diam_abs(maxDiameterY, (mpfi_ptr) x);
        mpfi_diam_abs(maxDiameterZ, (mpfi_ptr) x);
        
        //get max values between each of the mfpr widths
        mpfr_max(intermediateMax1, maxDiameterX, maxDiameterY, MPFR_RNDN);
        mpfr_max(intermediateMax2, intermediateMax1, maxDiameterZ, MPFR_RNDN);

        //compare mpfr max with the compare double value and return. if intermediateMax2 compare == 1, if they are equal, 
        //compare == 0, else compare == -1 so we need to correct for this in an if-else statement. 
        bool truth = mpfr_cmp_d(intermediateMax2, compare);

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

//function that checks if 0 is defined for both f and g
bool hasZero(mpfi_t x, mpfi_t y, mpfi_t z) {

    // initializing intervals to store f and g, and then computing f and g with the given x,y,z intervals. 
    mpfi_t resultOfF;
    mpfi_init((mpfi_ptr) resultOfF);

    mpfi_t resultOfg;
    mpfi_init((mpfi_ptr) resultOfg);

    f(resultOfF, x, y, z);

    g(resultOfg, x, y, z);

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
    intervalClass initialInterval = intervalClass(-1, 1, -1, 1, -1, 1);


    //define error
    double maxError = 1.0;




    //initialize queue of all intervals 
    std::queue<intervalClass> intervalQueue;
    //update with the first box
    intervalQueue.push(initialInterval);

    //initialize empty answer vector of intervals to contain where f and g are both equal to 0
    std::vector<intervalClass> answerIntervals;



    mpfr_t maxDiameterX;
    mpfr_t maxDiameterY;
    mpfr_t maxDiameterZ;

    mpfr_t intermediateMax1;
    mpfr_t intermediateMax2;

    mpfr_init(maxDiameterX);
    mpfr_init(maxDiameterY);
    mpfr_init(maxDiameterZ);

    mpfr_init(intermediateMax1);
    mpfr_init(intermediateMax2);

    mpfr_max((mpfr_ptr)intermediateMax1, maxDiameterX, maxDiameterY, MPFR_RNDN);































    // std::queue<student> q;
    // for(int i = 0; i< 5; i++){
    //     student newStudent = student(i);
    //     q.push(newStudent);
    // }
    // for(int i = 0; i< 5; i++){
    //     student newStudent = q.front();
    //     q.pop();
    //     std::cout << newStudent.score << &newStudent<< std::endl;
    // }

    // std::queue<int> queue;
    // for(int i = 0; i< 5; i++){
    //     queue.push(i);
    // }
    // for(int i = 0; i< 5; i++){
    //     int x = queue.front();
    //     queue.pop();
    //     std::cout << x << &x<< std::endl;
    // }

    // student arr[2];
    // for(int i = 0; i < 2; i++){
    //     arr[i] = student(i);
    // }
    // for(int i = 0; i < 2; i++){
    //     std::cout<< arr[i].score;
    // }
    

    //testing whether intervals have f and g

    // mpfi_t x;
    // // Initialize the interval variable
    // mpfi_init((mpfi_ptr) x);
    // // Set the interval to [-1, 1]
    // mpfi_interv_d( (mpfi_ptr) x, 0, 1);

    // mpfi_t y;
    // // Initialize the interval variable
    // mpfi_init((mpfi_ptr) y);
    // // Set the interval to [-2, 1]
    // mpfi_interv_d( (mpfi_ptr) y, 0, 1);

    // mpfi_t z;
    // // Initialize the interval variable
    // mpfi_init((mpfi_ptr) z);
    // // Set the interval to [-1, 1]
    // mpfi_interv_d( (mpfi_ptr) z, 0, 1);


    // //first initialize arrays
    // mpfi_t newX[2];
    // mpfi_t newY[2];
    // mpfi_t newZ[2];

    // //initializing mpfi_t arrays before placing values into them
    // for(int i = 0; i<2; i++){
    //     mpfi_init(newX[i]);
    //     mpfi_init(newY[i]);
    //     mpfi_init(newZ[i]);
    // }

    // //splitting box into 8 now boxes, so splitting x,y,z into 2 each so we can get 2 * 2 * 2 = 8 different combinations of boxes
    // mpfi_bisect((mpfi_ptr)newX[0], (mpfi_ptr) newX[1], x);
    // mpfi_bisect((mpfi_ptr)newY[0], (mpfi_ptr) newY[1], y);
    // mpfi_bisect((mpfi_ptr)newZ[0], (mpfi_ptr) newZ[1], z);


    // //triple for loop to queue all new elements. 
    // for(int i = 0; i <2 ; i++){
    //     for(int j = 0; j < 2; j++){
    //         for(int k = 0; k <2; k++){
    //             std::cout <<" intervals, X, Y, Z: " << std::endl;
    //             mpfi_out_str(stdout, 10, 0, (mpfi_ptr) newX[i]);
    //             std::cout<<std::endl;
    //             mpfi_out_str(stdout, 10, 0, (mpfi_ptr) newY[j]);
    //             std::cout<<std::endl;
    //             mpfi_out_str(stdout, 10, 0, (mpfi_ptr) newZ[k]);
    //             std::cout<<std::endl;
    //             std::cout << "Does the box contain it? "<<hasZero(newX[i], newY[j], newZ[k]) << std::endl;
    //             std::cout<<std::endl;
    //             std::cout<<std::endl;
    //             std::cout<<std::endl;
    //         }
    //     }
    // }

    // //clearing out variables that aren't needed
    // mpfr_t maxDiameterX;


    // mpfr_init(maxDiameterX);
    
    // // Set the variable to 0.25 using a double
    // //mpfr_set_d(maxDiameterX, 0.25, MPFR_RNDN);

    // mpfi_diam_abs(maxDiameterX, x);

    // bool truth = mpfr_cmp_d(maxDiameterX, 0.5);

    // std::cout<<std::endl;
    // std::cout << "Diameter of box? "<< std::endl;

    
    // std::cout << " Leq 0.25: "<< truth << std::endl;
    // std::cout<<std::endl;
    // std::cout<<std::endl;
    // //mpfi_out_str(stdout, 10, 0, (mpfi_ptr) maxDiameterX);
    // mpfr_out_str(stdout, 10, 0, maxDiameterX, MPFR_RNDN);
    // std::cout<<std::endl;

    // mpfr_clear(maxDiameterX);

    // for(int i = 0; i<2; i++){
    //     mpfi_clear(newX[i]);
    //     mpfi_clear(newY[i]);
    //     mpfi_clear(newZ[i]);
    // }

    // mpfi_clear((mpfi_ptr) x);
    // mpfi_clear((mpfi_ptr) y);
    // mpfi_clear((mpfi_ptr) z);



}