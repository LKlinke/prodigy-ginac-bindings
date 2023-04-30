#ifndef PRODIGY_ASSUMPTIONS_H
#define PRODIGY_ASSUMPTIONS_H

#include <string>
#include <memory>

namespace prodigy {
    class Assumption {
        public:
        virtual bool implies(std::shared_ptr<Assumption>) = 0;
    };

    class Divisible : public Assumption {
        private:
        std::string lhs;
        uint32_t rhs;
        public:
        Divisible(const std::string& lhs, const uint32_t& rhs);
        bool implies(std::shared_ptr<Assumption>) override;
        std::string getLhs() const;
        uint32_t getRhs() const;
    };
}

#endif

