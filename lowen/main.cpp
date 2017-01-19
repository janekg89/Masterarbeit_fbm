//
// Created by janekg89 on 11.07.16.
//

#include "_generatefracincrements.h"

int main(){
    Increments1 inc(100, 1);
    inc.generateIncrements1(2.0,0.5,0.7);
    boost::multi_array<double,3> data = inc.increments;
    for (auto i=0; i < data.shape()[2]; ++i) {
        std::cout << data[0][0][i] << std::endl;

    }
    return 0;
};