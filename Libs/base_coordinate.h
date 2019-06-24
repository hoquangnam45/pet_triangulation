#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <string.h>

class baseCoordinate{
    public:
        // constructor
        baseCoordinate(double x, double y, double z);
        baseCoordinate();
        virtual ~baseCoordinate();
        // set
        void setCoordinate(int x, int y, int z);
        // get
        double getX();
        double getY();
        double getZ();
        // copy
        baseCoordinate& operator=(const baseCoordinate& b);
        // read;
        friend std::istream & operator >> (std::istream &in,  baseCoordinate &c);
        // output
        friend std::ostream & operator << (std::ostream &out,  const baseCoordinate &c);
    protected:
        double x, y, z;
    private:
        void parseBuf(const char*);
        void printOut();
};