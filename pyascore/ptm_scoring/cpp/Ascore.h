#ifndef ASCORE_H
#define ASCORE_H

namespace ptmscoring {
    class Ascore {
        public:
            float window_size;
            Ascore();
            Ascore(float);
            ~Ascore();
            void getWindowSize();
    };
}

#endif
