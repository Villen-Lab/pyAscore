#include <iostream>
#include "Ascore.h"

namespace ptmscoring {

    Ascore::Ascore () {}

    Ascore::Ascore (float window_size) {
        this->window_size = window_size;
    }

    Ascore::~Ascore () {}

    void Ascore::getWindowSize () {
        std::cout << window_size << std::endl;
    }
}
