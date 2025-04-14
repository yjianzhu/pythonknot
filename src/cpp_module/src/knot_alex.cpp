#include "knot_alex_table.h"
#include <fstream>

GiNaC::symbol t("t");
GiNaC::lst syms={t};
std::map<GiNaC::ex,std::string,GiNaC::ex_is_less>  alexander_polynomial = 
{
    {GiNaC::ex("1",syms), "1"},
    {GiNaC::ex("-1+t-t^2",syms), "3_1"},
    {GiNaC::ex("-1+3*t-t^2",syms), "4_1"},
    {GiNaC::ex("1-t+t^2-t^3+t^4",syms), "5_1"},
    {GiNaC::ex("2-3*t+2*t^2",syms), "5_2"},
    {GiNaC::ex("-2+5*t-2*t^2",syms), "6_1"},
    {GiNaC::ex("1-3*t+3*t^2-3*t^3+t^4",syms), "6_2"},
    {GiNaC::ex("1-3*t+5*t^2-3*t^3+t^4",syms), "6_3"},
    {GiNaC::ex("1-t+t^2-t^3+t^4-t^5+t^6",syms), "7_1"},
    {GiNaC::ex("3-5*t+3*t^2",syms), "7_2"},
    {GiNaC::ex("2-3*t+3*t^2-3*t^3+2*t^4",syms), "7_3"},
    {GiNaC::ex("4-7*t+4*t^2",syms), "7_4"},
    {GiNaC::ex("2-4*t+5*t^2-4*t^3+2*t^4",syms), "7_5"},
    {GiNaC::ex("1-5*t+7*t^2-5*t^3+t^4",syms), "7_6"},
    {GiNaC::ex("1-5*t+9*t^2-5*t^3+t^4",syms), "7_7"},
    {GiNaC::ex("3-7*t+3*t^2",syms), "8_1"},
    {GiNaC::ex("1-3*t+3*t^2-3*t^3+3*t^4-3*t^5+t^6",syms), "8_2"},
    {GiNaC::ex("4-9*t+4*t^2",syms), "8_3"},
    {GiNaC::ex("2-5*t+5*t^2-5*t^3+2*t^4",syms), "8_4"},
    {GiNaC::ex("1-3*t+4*t^2-5*t^3+4*t^4-3*t^5+t^6",syms), "8_5"},
    {GiNaC::ex("2-6*t+7*t^2-6*t^3+2*t^4",syms), "8_6"},
    {GiNaC::ex("1-3*t+5*t^2-5*t^3+5*t^4-3*t^5+t^6",syms), "8_7"},
    {GiNaC::ex("2-6*t+9*t^2-6*t^3+2*t^4",syms), "8_8"},
    {GiNaC::ex("1-3*t+5*t^2-7*t^3+5*t^4-3*t^5+t^6",syms), "8_9"},
    {GiNaC::ex("1-3*t+6*t^2-7*t^3+6*t^4-3*t^5+t^6",syms), "8_10"},
    {GiNaC::ex("2-7*t+9*t^2-7*t^3+2*t^4",syms), "8_11"},
    {GiNaC::ex("1-7*t+13*t^2-7*t^3+t^4",syms), "8_12"},
    {GiNaC::ex("2-7*t+11*t^2-7*t^3+2*t^4",syms), "8_13"},
    {GiNaC::ex("2-8*t+11*t^2-8*t^3+2*t^4",syms), "8_14"},
    {GiNaC::ex("3-8*t+11*t^2-8*t^3+3*t^4",syms), "8_15"},
    {GiNaC::ex("1-4*t+8*t^2-9*t^3+8*t^4-4*t^5+t^6",syms), "8_16"},
    {GiNaC::ex("1-4*t+8*t^2-11*t^3+8*t^4-4*t^5+t^6",syms), "8_17"},
    {GiNaC::ex("1-5*t+10*t^2-13*t^3+10*t^4-5*t^5+t^6",syms), "8_18"},
    {GiNaC::ex("1-t+t^3-t^5+t^6",syms), "8_19"},
    {GiNaC::ex("1-2*t+3*t^2-2*t^3+t^4",syms), "8_20"},
    {GiNaC::ex("1-4*t+5*t^2-4*t^3+t^4",syms), "8_21"}
};

void get_alexander(std::fstream &read)
{   
    // clear the map, and the insert // or 先检查再insert
    alexander_polynomial.clear();
    
    std::string knot_type,alex_poly;
    while ( !read.eof())
    {
        read>>knot_type>>alex_poly;
        GiNaC::ex alex_poly_ex(alex_poly,syms);
        //cout<<alex_poly_ex<<'\t'<<knot_type<<endl;
        
        alexander_polynomial.insert(std::pair<GiNaC::ex,std::string>(alex_poly_ex,knot_type));
        
        // if(alexander_polynomial.count(alex_poly_ex)==0)
        // {
        //     alexander_polynomial.insert(std::pair<GiNaC::ex,std::string>(alex_poly_ex,knot_type));
        // }
    }
}

void print_alexander_map()
{
    std::cout<<alexander_polynomial.size()<<std::endl;
    for(auto t=alexander_polynomial.begin();t!=alexander_polynomial.end();t++){

    std::cout<<t->first<<'\t'<<t->second<<std::endl;

    }
}