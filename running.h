#ifndef RUNNING_H
#define RUNNING_H

class Running
{
public:

    char field_type;
    int maxit, systemsize;  // Should I be allowed to vary J during a set of runs?
    double tolerance, J;
    bool armadillobool, sectorbool;

    Running(); // Do I even need this class? No, I don't think I do...
};

#endif // RUNNING_H
