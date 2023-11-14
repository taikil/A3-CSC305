#include <iostream>
#include <cmath>

using namespace std;

int main()
{
    double sales = 95000;
    double sTax = 0.04;
    double cTax = 0.02;
    cout << "Sales: $" << sales << endl
         << "State Tax: $" << sales * sTax << endl
         << "County Tax: $" << sales * cTax << "\n";
    return 0;
}