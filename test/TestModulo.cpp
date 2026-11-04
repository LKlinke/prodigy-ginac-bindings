// Sections of rational gfs modulo q (arithmetic progressions, issue #10): they must be rational functions with
// rational coefficients, without roots of unity, and they must have the same power series as the sections of the gf.

#include <DistributionFactory.h>
#include "gtest/gtest.h"

using namespace prodigy;
using namespace GiNaC;

namespace {

    std::string str(const Dist &d) {
        std::stringstream s;
        s << d;
        return s.str();
    }

    // [x^k] of the power series of e in x, computed with GiNaC's series, independently of the multisection
    ex seriesCoeff(const ex &e, const symbol &x, int k) {
        return series_to_poly(e.series(x == 0, k + 1)).expand().coeff(x, k).normal();
    }

    void expectSections(const std::string &gf, const std::vector<std::string> &params, int q, bool rational = true,
                        int order = 13) {
        Dist dist(gf, params);
        std::vector<Dist> sections = dist.arithmeticProgression("x", std::to_string(q));
        ASSERT_EQ(sections.size(), static_cast<std::size_t>(q));

        // the parser of the test uses its own symbols, so the comparison does not depend on Dist's symbol table
        symtab table;
        symbol x("x");
        table["x"] = x;
        parser reader(table);
        const ex f = reader(gf);

        Dist sum("0", params);
        for (int r = 0; r < q; r++) {
            const std::string s = str(sections[r]);
            if (rational) {
                EXPECT_EQ(s.find("Pi"), std::string::npos) << gf << ", q = " << q << ", r = " << r << ": " << s;
                EXPECT_EQ(s.find('I'), std::string::npos) << gf << ", q = " << q << ", r = " << r << ": " << s;
            }
            const ex section = reader(s);
            for (int k = 0; k < order; k++) {
                const ex expected = k % q == r ? seriesCoeff(f, x, k) : ex(0);
                const ex diff = (seriesCoeff(section, x, k) - expected).normal();
                if (rational) {
                    EXPECT_TRUE(diff.is_zero()) << gf << ", q = " << q << ", r = " << r << ", k = " << k << ": " << s;
                } else {
                    // GiNaC cannot decide identities of roots of unity (issue #10), so this test of the fallback
                    // compares the coefficients numerically
                    EXPECT_TRUE(abs(diff).evalf() < numeric("1e-12"))
                        << gf << ", q = " << q << ", r = " << r << ", k = " << k << ": " << s;
                }
            }
            sum = sum + sections[r];
        }
        if (rational) {
            EXPECT_EQ(sum, dist) << gf << ", q = " << q;
        }
    }
}

TEST(Modulo, issue10Examples) {
    std::vector<Dist> three = Dist("1/2/(1-1/2*x)").arithmeticProgression("x", "3");
    EXPECT_EQ(three[0], Dist("1/(2-x^3/4)"));
    EXPECT_EQ(three[1], Dist("x/2/(2-x^3/4)"));
    EXPECT_EQ(three[2], Dist("x^2/4/(2-x^3/4)"));
    std::vector<Dist> five = Dist("1/2/(1-1/2*x)").arithmeticProgression("x", "5");
    EXPECT_EQ(five[0], Dist("1/(2-x^5/16)"));
    EXPECT_EQ(five[3], Dist("x^3/8/(2-x^5/16)"));
    std::vector<Dist> two = Dist("1/2/(1-1/2*x)").arithmeticProgression("x", "2");
    EXPECT_EQ(two[0], Dist("1/(2-x^2/2)"));
}

TEST(Modulo, univariate) {
    for (int q : {1, 2, 3, 4, 5, 6, 7}) {
        expectSections("1/2/(1-1/2*x)", {}, q);
        expectSections("(1+x^2)/((2-x)*(3-x))", {}, q);
        expectSections("1/(1-x/3-x^2/3)", {}, q);
        expectSections("(1/4+x/4)^3*(1+x^7)", {}, q);
        expectSections("x^4/(2-x)^3", {}, q);
        expectSections("1/9*(x^5+2)/(1-x/3)^2/(1+x^2/5)", {}, q);
    }
}

TEST(Modulo, multivariate) {
    for (int q : {2, 3, 5, 6}) {
        expectSections("1/(3-x-y)", {}, q);
        expectSections("y*(1+x*y)/((2-x*y)*(4-x-y^2))", {}, q);
        expectSections("1/2/(1-1/2*x)*z^2*(1/2+1/2*y)", {}, q);
    }
}

TEST(Modulo, parameters) {
    for (int q : {2, 3, 5, 6}) {
        expectSections("p/(1-(1-p)*x)", {"p"}, q);
        expectSections("(p*x+1-p)^4*y", {"p"}, q);
        expectSections("(1-p)*(1+a*x)/(1-p*x^2*y)", {"p", "a"}, q);
        expectSections("exp(-l)/(2-x)", {"l"}, q);
        expectSections("sqrt(p)/(2-sqrt(p)*x)", {"p"}, q);
    }
}

TEST(Modulo, notRationalInX) {
    // gfs that are not rational in x keep the roots of unity, but their sections stay exact
    for (int q : {2, 3}) {
        expectSections("exp(x-1)", {}, q, false, 8);
        expectSections("exp(x-1)/(2-x)", {}, q, false, 8);
    }
}

TEST(Modulo, updateModulo) {
    // y := x % 3 on a geometric x
    Dist res = Dist("1/2/(1-1/2*x)").updateModulo("y", "x", "3", {"x", "y"}, false, std::nullopt);
    EXPECT_EQ(res, Dist("(1+x*y/2+x^2*y^2/4)/(2-x^3/4)")) << res;
    std::string s = str(res);
    EXPECT_EQ(s.find("exp"), std::string::npos) << s;
}
