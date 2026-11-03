#include "Assumptions.h"

namespace prodigy {
    Divisible::Divisible(const std::string& lhs, const uint32_t& rhs) {
        this->lhs = lhs;
        this->rhs = rhs;
    }

    std::string Divisible::getLhs() const {
        return this->lhs;
    }

    uint32_t Divisible::getRhs() const {
        return this->rhs;
    }

    bool Divisible::implies(std::shared_ptr<Assumption> other) {
        auto as_divisible = dynamic_cast<Divisible*>(other.get());
        if(as_divisible == nullptr) {
            return false;
        }
        return as_divisible->lhs == this->lhs && this->rhs % as_divisible->rhs == 0;
    }
}