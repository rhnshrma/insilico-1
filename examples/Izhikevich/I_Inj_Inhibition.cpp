#include<iostream>

#include<fstream>
using namespace std;
int main() 
{
	double a,b;// a is time and b is current
	ofstream out_data("Inhibition_induced_spiking_data.isf");
	out_data <<"time,0"<<endl;
	for (double a = 100,b = 0; a <=801;a = a + 0.01){
		b = -50 ;
		out_data <<a<<","<<b<<endl;
		}
	return 0;
}

