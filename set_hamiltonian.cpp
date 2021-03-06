#include "set_hamiltonian.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <typeinfo>

using namespace std;

Set_Hamiltonian::Set_Hamiltonian()
{

}

Set_Hamiltonian::Set_Hamiltonian(int systemsize, double J, vector<double> hs, bool armadillobool, bool sectorbool)
{
    this->systemsize = systemsize;                // Is this notation really neccessary? Look it up
    this->no_of_states = pow(2, systemsize);
    this->J=J;
    this->hs=hs;
    this->armadillobool = armadillobool;
    this->sectorbool = sectorbool;
    palhuse = true;                               // Default setting
    //cout << "In constructor" << endl;
    // Should do things more automatically, probably... Or?
}

void Set_Hamiltonian::donotconservespin()         // I don't need this yet. It might be useful later
{
    palhuse = false;
}

void Set_Hamiltonian::give_hxs(vector<double> hxs_in)  // For accessing the Set_Hamiltonian class directly (for testing and such)
{
    this->hxs = hxs_in;
    cout << "hxs given" << endl;
    /*
    cout << "Printing hxs: " << endl;
    for(int i=0; i<systemsize; i++)
    {
        cout << "hxs_in[" << i << "] = " << hxs_in[i] << endl;
        cout << "hxs[" << i << "] = " << hxs[i] << endl;
    }
    */

}

void Set_Hamiltonian::create_armadillo_matrix()
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    cout << "In create_armadillo_matrix" << endl;
    armaH = arma::mat(no_of_states, no_of_states);
    cout << "Matrix initialized" << endl;
    for(int i=0; i<no_of_states; i++)
    {
        for(int j=0; j<no_of_states; j++)    armaH(i,j)= 0.0;
    }
    //cout << "Arma matrix set. Max index no = " << no_of_states-1 << endl;
}

void Set_Hamiltonian::create_armadillo_matrix(int size)
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    armaH = arma::mat(size, size);
    for(int i=0; i<size; i++)        // This does not seem to be neccessary all the time, but it is safer.
    {
        for(int j=0; j<size; j++)    armaH(i,j)= 0.0;
    }
    //cout << "Arma matrix set. Max index no = " << size-1 << endl;
}


void Set_Hamiltonian::create_dense_Eigen_matrix()
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    //const bool TRACE = false;
    //if(TRACE)    cout << "I am going to set the size of an eigenmatrix" << endl;
    eigenH = Eigen::MatrixXd(no_of_states, no_of_states);
    for(int i=0; i<no_of_states; i++)
    {
        for(int j=0; j<no_of_states; j++)    eigenH(i,j)= 0.0;
    }
    //cout << "Eigen matrix set. Max index no. = " << no_of_states-1 << endl;
}

void Set_Hamiltonian::create_dense_Eigen_matrix(int size)
{   // Creates an empty matrix, ready to use, in those cases where it is not too big.
    //const bool TRACE = false;
    //if(TRACE)    cout << "I am going to set the size of an eigenmatrix" << endl;
    eigenH = Eigen::MatrixXd(size, size);
    for(int i=0; i<size; i++)
    {
        for(int j=0; j<size; j++)    eigenH(i,j)= 0.0;
    }
    //cout << "Eigen matrix set. Max index no. = " << size-1 << endl;
}





//---------------------------------------BASIC SPIN OPERATIONS--------------------------------------------//

int Set_Hamiltonian::give_spin(int i, int a)
{
    return ((a&(1<<i))>>i);
}

int Set_Hamiltonian::set_spin_up(int i, int a)
{
    return (a |(1<<i));
}

int Set_Hamiltonian::set_spin_down(int i, int a)
{
    return ~((~a)|(1<<i));
}

int Set_Hamiltonian::flip_spin(int i, int a)      // Used only by sxi
{
    return (a ^(1<<i));
}


//--------------------------------------------SPIN OPERATORS----------------------------------------------//
//-------------------------------------------------/Sz/---------------------------------------------------//


double Set_Hamiltonian::szi(int i, int a)
{
    if(give_spin(i,a)==0)    return -0.5;
    else                     return 0.5;
}

double Set_Hamiltonian::sxi(int i, int a)
{   // So the idea here is: S^x = 0.5*(S^+ + S^-). When acting on a state,
    // one of these terms will flip the i'th spin, the other will return
    //nothing. So this is the same as flipping the spin in the i-position
    return 0.5*flip_spin(i,a);
}

double Set_Hamiltonian::szip1szi(int i, int a)
{
    double firstspin = szi(i, a);
    double secondspin = 0;
    if(i==(systemsize-1))    secondspin = szi(0,a);
    else                     secondspin = szi((i+1),a);
    return firstspin*secondspin;
}


//-------------------------------------------/S+ and S-/--------------------------------------------------//


int Set_Hamiltonian::upip1downi(int i, int a)
{
    int a2 = 0;
    if(i==(systemsize-1))            a2 = set_spin_up(0, a);         // Periodic boundary conditions
    else if(i<(systemsize-1))        a2 = set_spin_up((i+1), a);     //
    else                             cout << "Check your indices, woman!!" << endl;
    return set_spin_down(i, a2);
}

int Set_Hamiltonian::downip1upi(int i, int a)
{
    int a2 = 0;
    if(i==(systemsize-1))            a2 = set_spin_down(0, a);       // Periodic boundary conditions
    else if(i<(systemsize-1))        a2 = set_spin_down((i+1), a);
    else                             cout << "Check your indices, woman!!" << endl;
    return set_spin_up(i, a2);
}


//-----------------------------------------STATE-SPECIFIC FUNCTIONS---------------------------------------//
// Both of these functions work. It might be excessive to have number_of_down(a,systemsize) as a separate
// funtction.
int Set_Hamiltonian::number_of_up(int a)
{
    int no_of_up = 0;
    for(int i=0; i<systemsize; i++)    no_of_up += give_spin(i, a);
    return no_of_up;
}

int Set_Hamiltonian::number_of_down(int a)
{
    return (systemsize - number_of_up(a));
}

// function for s0 sector

void Set_Hamiltonian::checktestupdown(int j,int i)
{
    int helpnumber = 0;


    if(j!=(systemsize-1))
    {
        testupip1downi = (give_spin(j,i)==1) && ( (give_spin((j+1),i))==0);
        testdownip1upi = (give_spin(j,i)==0) && ( (give_spin((j+1),i))==1);
    }
    else
    {
        testupip1downi = (give_spin(j,i)==1) && ( (give_spin(helpnumber,i))==0);
        testdownip1upi = (give_spin(j,i)==0) && ( (give_spin(helpnumber,i))==1);
    }
}

/*
void Set_Hamiltonian::set_mysector(int mysector)
{
    if(sectorbool==false)  // Just in case I am that stupid.
    {
        sectorbool = true;
        sectorboolchanged = true;
        cout << "NB! Sectorbool changed to true!" << endl;
    }
    this->mysector = mysector;
    find_sector_dense();
}
*/

void Set_Hamiltonian::set_mysector(int mysector)
{   // Should keep some information on what states are in the list. Well, that is all in sectorlist.
    // And that is not destroyed. Retrieve it in main, perhaps.

    // Initializing
    no_of_hits = 0;
    sectorlist = vector<int>(no_of_states);
    this->mysector = mysector;

    // A warning will probably be nice...
    if(sectorbool==false)  // Just in case I am that stupid.
    {
        sectorbool = true;
        sectorboolchanged = true;
        cout << "NB! Sectorbool changed to true!" << endl;
    }

    // Actually finding the sector
    for(int state=0; state<no_of_states; state++)
    {
        if(number_of_up(state)==mysector)
        {
            sectorlist[no_of_hits] = state;
            no_of_hits ++;
        } // End if
    } // End for
    if(no_of_hits > 0)                   trim_sectorlist();
    else
    {
        sectorlist = vector<int>(1);
        sectorlist[0] = 0;
        no_of_hits = 1;
    }
}

void Set_Hamiltonian::trim_sectorlist()
{   // A function that cuts sectorlist off as it detects a zero
    vector<int> sectorlist_short = vector<int>(no_of_hits);
    for(int i=0; i<no_of_hits; i++)     sectorlist_short[i] = sectorlist[i];
        sectorlist = sectorlist_short;
}

//--------------------------------------------//THE SYSTEM//----------------------------------------------//

//-------------------------------------------/HAMILTONIANS/--------------------------------------------//

void Set_Hamiltonian::set_elements(int i, int b)    // Do I really need this function? --- Only if I am going to have some non-PalHuse Hamiltonian with spin conservation.
{
    if(palhuse==true)    palhuse_set_elements(i, b);
}

void Set_Hamiltonian::palhuse_set_elements(int i, int b)
{
    cout << "In palhuse_set_elements" << endl;
    int index1 = 0;
    int index2 = 0;
    if(sectorbool==false)
    {
        index1 = no_of_states - (b+1);  // Is this really neccessary? The matrix is symmetric, so the off-diag part will be correct either way.
        index2 = no_of_states - (i+1);
    }
    else
    {   // Must correct this later. Do not want an upside-down matrix. Or? Not a big problem if I retrieve
        // sectorlist. That holds for all systems with the same no. of. particles.
        index1 = i;
        index2 = b;
        // Blimey. This redistribution of states is trickier than I thought.
    }


    double element = 0;
    if(systemsize == 2)        element = J;  // Special case. BCs and small
    else                       element = 0.5*J;

    if(armadillobool == true)
    {
        armaH(index1,index2) = element;
        armaH(index2,index1) = element;
    }
    else
    {
        eigenH(index1, index2) = element;
        eigenH(index2, index1) = element;

    }
}

void Set_Hamiltonian::set_element_spinnnotconserved(int i, int j)
{
    cout << "In set_element_spinnotconserved" << endl;
    int b = flip_spin(j,i);
    if(armadillobool)        armaH(b,i) = 0.25*hxs[j];
    else                     eigenH(b,i) = 0.25*hxs[j];

}


//----------------------------------------SECTOR HAMILTONIANS---------------------------------------------//

void Set_Hamiltonian::palhuse_interacting_sectorHamiltonian()
{
    const bool TRACE = false;
    if(TRACE)    cout << "At least I am in palhuse_interacting_sectorHamiltonian" << endl;

    matrixsize = no_of_hits;

    if(armadillobool)        create_armadillo_matrix(no_of_hits);
    else                     create_dense_Eigen_matrix(no_of_hits);

    int a = 0;
    int b = 0;
    for(int i=0; i<no_of_hits; i++)
    {
        a = sectorlist[i];
        for(int j=0; j<systemsize; j++)
        {   // Should I include a while loop, or something?
            checktestupdown(j,a);
            if(testupip1downi==true)
            {
                b = upip1downi(j, a);
                for(int k = 0; k<no_of_hits; k++)
                {
                    if(b==sectorlist[k])  set_elements(i, k);
                }
            }
            if(testdownip1upi==true)
            {
                b = downip1upi(j, a);
                for(int k = 0; k<no_of_hits; k++)
                {
                    if(b==sectorlist[k])  set_elements(i, k);
                }
            }
        } // End for j
    } // End for i
    if(TRACE)    cout << "Exiting palhuse_interacting_sectorHamiltonian" << endl;
}

void Set_Hamiltonian::palhuse_diagonal_sectorHamiltonian()
{
    const bool TRACE = false;
    // Do something like index = no_of_states - i that works for this one.
    // For now, the matrix is upside down, but we have gotten a list of its entries: sectorlist.
    double element = 0;
    //double index;
    for(int i=0; i<no_of_hits; i++)
    {
        element = 0;
        //index = no_of_hits-1-i;  // To order the matrix from highest to lowest number in the binary representation
        int a = sectorlist[i];
        for(int j=0; j<systemsize; j++)  element += hs[j]*szi(j, a) + J*szip1szi(j,a);
        if(armadillobool == true)                            armaH(i,i) = element;
        else if(armadillobool == false)                      eigenH(i,i) = element;
    } // End for-loop over i
    if(TRACE)    cout << "Exiting palhuse_diagonal_sectorHamiltonian" << endl;
} // End function palhuse_random_sectorHamiltonian_dense

//-------------------------------------------TOTAL HAMILTONIAN--------------------------------------------//

void Set_Hamiltonian::palhuse_interacting_totalHamiltonian()
{
    if(armadillobool)                                               create_armadillo_matrix();
    else                                                            create_dense_Eigen_matrix();
    sectorbool=false; // Is this neccessary?

    matrixsize = no_of_states;

    int b = 0;
    for(int i=0; i<no_of_states; i++)
    {   // i is our state
        for(int j=0; j<systemsize; j++)
        {
            checktestupdown(j,i);
            if(testupip1downi==true)
            {
                b = upip1downi(j, i);
                set_elements(i,b);
            }
            if(testdownip1upi==true) // else if here would maybe save some time (not much, though)
            {
                b = downip1upi(j, i);
                set_elements(i,b);
            }
        } // End for j
    } // End for i
} // End function


void Set_Hamiltonian::palhuselike_spinnotconserved_interacting_totalHamiltonian()
{
    //cout << "In palhuselike_spinnotconserved_interactiong_totalHamiltonian" << endl;
    if(armadillobool)                                               create_armadillo_matrix();
    else                                                            create_dense_Eigen_matrix();
    sectorbool=false; // Is this neccessary?

    matrixsize = no_of_states;

    int b = 0;
    for(int i=0; i<no_of_states; i++)
    {   // i is our state
        for(int j=0; j<systemsize; j++)
        {
            cout << "in loop" << endl;
            // Setting the elements that we did in Pal & Huse
            checktestupdown(j,i);
            if(testupip1downi==true)
            {
                cout << "S+S-" << endl;
                b = upip1downi(j, i);
                set_elements(i,b);
            }
            if(testdownip1upi==true) // else if here would maybe save some time (not much, though)
            {
                cout << "S-S+" << endl;
                b = downip1upi(j, i);
                set_elements(i,b);
            }
            if(armadillobool)
            {
                for(int i=0; i<no_of_states; i++)        // This does not seem to be neccessary all the time, but it is safer.
                {
                    for(int j=0; j<no_of_states; j++)    cout << armaH(i,j) << " ";
                    cout << endl;
                }
            }
            else
            {
                for(int i=0; i<no_of_states; i++)        // This does not seem to be neccessary all the time, but it is safer.
                {
                    for(int j=0; j<no_of_states; j++)    cout << eigenH(i,j) << " ";
                    cout << endl;
                }
            }
            //cout << "If-tests ended for now. Setting the sxi-generated elements: " << endl;
            //set_element_spinnnotconserved(i,j);
        } // End for j
    } // End for i
}



void Set_Hamiltonian::palhuse_diagonal_totalHamiltonian()
{
    cout << "In palhuse_diagonal_totalHamiltonian()" << endl;
    int index = 0;
    double element = 0;
    cout << "Just testing to see if, for some stupid reason, I can't access the elements of a vector:" << endl;
    cout << "Here goes nothing: hs[0] = " << hs[0] << endl;
    cout << "You're here, so the fault must be somewhere else" << endl;

    cout << "Before setting the diagonal elements" << endl;
    for(int i=0; i<no_of_states; i++)   // Loop over every matrix element
    {
        element = 0;
        // The 2-particle case is special again...
        for(int j=0; j<systemsize; j++)  element += hs[j]*szi(j, i) + J*szip1szi(j,i);   // Loop over the contributions from each spin site.
        index = no_of_states - (i+1);
        //index = i;
        if(armadillobool==true)                            armaH(index,index) = element;
        else if(armadillobool==false)                      eigenH(index,index) = element;
    } // End for-loop over i
    cout << "After setting the diagonal elements" << endl;
    if(armadillobool==false)                                                  eigenH.cast<double>(); // Is this really neccessary?
} // End function palhuse_random_sectorHamiltonian_dense



// Do I need something like this for sectors, maybe?
/*

void Set_Hamiltonian::flip_diagonal_sectorHamiltonian()
{
    // See if I have to do something to flip the Hamiltonian. Could probably do this without calling its own function.

}

*/
