//
// Created by Tobias Winkler on 01.11.21.
//

#include "gtest/gtest.h"
#include "Dist.h"
#include "DistributionFactory.h"

using namespace std;
using namespace prodigy;

TEST(EquivalenceCheck, geometric) {
    Dist l0, l1, l2, l3, l4, l5, l6, l7;

    /*
     * given loop
     *
     * while(x = 1) {
     *     { x := 0 } [1/2] { c := c+1 }
     * }
     */

    /*
     * invariant
     *
     0:
     1: if(x = 1) {
     2:     c := c + geometric(1/2)
     3:     x := 0
     4: }
     *
     */
    l0 = DistributionFactory::geometric("x", "p") * DistributionFactory::geometric("c", "q");
    l1 = l0.filterEq("x", "1");
    l2 = l1 * DistributionFactory::geometric("c", "1/2");
    l3 = l2.update("x", "0");
    l4 = l3 + (l0 - l1);
    Dist res1 = l4;

    /*
     * guard + body + invariant
     *
     0:
     1: if(x = 1) {
     2:     { x := 0 } [1/2] { c := c+1 }
     3:     if(x = 1) {
     4:         c := c + geometric(1/2)
     5:         x := 0
     6:     }
     7: }
     *
     */
    l0 = DistributionFactory::geometric("x", "p") * DistributionFactory::geometric("c", "q");
    l1 = l0.filterEq("x", "1");
    l2 = "1/2" * l1.update("x", "0") + "1/2" * l1.update("c", "c + 1");
    l3 = l2.filterEq("x", "1");
    l4 = l3 * DistributionFactory::geometric("c", "1/2");
    l5 = l4.update("x", "0");
    l6 = l5 + (l2 - l3);
    l7 = l6 + (l0 - l1);
    Dist res2 = l7;

    ASSERT_EQ(res1, res2);
}

TEST(EquivalenceCheck, geometricSum) {
    Dist::resetSymbolCache();
    Dist l0, l1, l2, l3, l4, l5, l6, l7, l8;

    /*
     * given loop
     *
     * while(x > 0) {
     *     { x := x-1 } [1/2] { c := c+1 }
     * }
     *
     */

    /*
     * invariant
     0:
     1: tmp := iid(geometric(1/2), x)
     2: c := c + tmp
     3: tmp := 0
     4: x := 0
     *
     */
    l0 = DistributionFactory::geometric("x", "p") * DistributionFactory::geometric("c", "q");
    // l0 = Dist {"1/(1 - p*x)", "p"} * Dist {"1/ (1 - q*c)", "q"};
    l1 = l0.updateIid("tmp", DistributionFactory::geometric("TMP", "1/2"), "x");
    l2 = l1.update("c", "c + tmp");
    l3 = l2.update("tmp", "0");
    l4 = l3.update("x", "0");
    Dist res1 = l4;

    /*
     * guard + body + invariant
     *
     0:
     1: if(x > 0) {
     2:     { x := x-1 } [1/2] { c := c+1 }
     3:     tmp := iid(geometric(1/2), x)
     4:     c := c + tmp
     5:     tmp := 0
     6:     x := 0
     7: }
     *
     *
     */
    l0 = DistributionFactory::geometric("x", "p") * DistributionFactory::geometric("c", "q");
    // l0 = Dist {"1/(1 - p*x)", "p"} * Dist {"1/ (1 - q*c)", "q"};
    l1 = l0.filterGreater("x", "0");
    l2 = "1/2" * l1.update("x", "x - 1") + "1/2" * l1.update("c", "c + 1");
    l3 = l2.updateIid("tmp", DistributionFactory::geometric("TMP", "1/2"), "x");
    l4 = l3.update("c", "c + tmp");
    l5 = l4.update("tmp", "0");
    l6 = l5.update("x", "0");
    l7 = l6 + (l0 - l1);
    Dist res2 = l7;

    ASSERT_EQ(res1, res2);
}

TEST(EquivalenceCheck, randomWalk) {
    Dist::resetSymbolCache();
    Dist l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16;

    /*
     * given loop
     *
     0:
     1: while (s > 0) {
     2:     { s:= s + 1 } [1/2] { s:= s - 1 }
     3:     c := c + 1
     5:     tmp := 0
     6: }
     *
     */

    /*
     * invariant
     *
     0:
     1: if(s > 0) {
     2:     tmp := iid(Dist{ (1 - sqrt(1-c^2)) / c }, s)
     3:     c := c + tmp
     4:     tmp := 0
     5:     s := 0
     6: }
     */
    l0 = Dist {"1/(1 - p*s)", "p"} *  Dist {"1/(1 - q*c)", "q"} * Dist {"1/(1 - r*tmp)", "r"};
    l1 = l0.filterGreater("s", "0");
    l2 = l1.updateIid("tmp", Dist {"(1 - (1-T^2)^(1/2)) / T"}, "s" );
    l3 = l2.update("c", "c + tmp");
    l4 = l3.update("tmp", "0");
    l5 = l4.update("s", "0");
    l6 = l5 + l0 - l1;
    Dist res1 = l6;

    /*
     * guard + body + invariant
     *
     0:
     1: if (s > 0) {
     2:     { s:= s+1 } [1/2] { s:= s-1 }
     3:     c := c + 1
     4:     tmp := 0
     5:     if(s > 0) {
     6:         tmp := iid(Dist{ (1 - sqrt(1-c^2)) / c }, s)
     7:         c := c + tmp
     8:         tmp := 0
     9:         s := 0
    10:     }
    11: }
     *
     */
    l0 = Dist {"1/(1 - p*s)", "p"} *  Dist {"1/(1 - q*c)", "q"} * Dist {"1/(1 - r*tmp)", "r"};
    l1 = l0.filterGreater("s", "0");
    l2 = "1/2" * l1.update("s", "s + 1") + "1/2" * l1.update("s", "s - 1");
    l3 = l2.update("c", "c + 1");
    l4 = l3.update("tmp", "0");
    l5 = l4.filterGreater("s", "0");
    l6 = l5.updateIid("tmp", Dist {"(1 - (1-T^2)^(1/2)) / T"}, "s" );
    l7 = l6.update("c", "c + tmp");
    l8 = l7.update("tmp", "0");
    l9 = l8.update("s", "0");
    l10 = l9 + l4 - l5;
    l11 = l10 + l0 - l1;
    Dist res2 = l11;

    ASSERT_EQ(res1, res2);
}

TEST(EquivalenceCheck, trivialIID){
    Dist::resetSymbolCache();
    Dist l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12;
    /*
    * given loop
    *
    0:
    1: while (0 < n) {
    2:   tmp := unif(1,6);
    3:   m := m + tmp;
    4:   tmp := 0;
    5:   n := n-1
    6: }
    *
    */

    /*
     * invariant
     *
     0:
     1: if (0 < n){
     2:   tmp := iid(unif(1,6), n)
     3:   m := m+tmp
     4:   n := 0
     5:   tmp := 0
     6:}
     */
    cout << "Invariant" << endl;
    l0 = Dist {"1/(1 - p*n)", "p"} *  Dist {"1/(1 - q*m)", "q"} * Dist {"1/(1 - r*t)", "r"};
    cout << "l0:\t" << l0 << endl;
    l1 = l0.filterGreater("n", "0");
    cout << "l1:\t" << l1 << endl;
    l2 = l1.updateIid("t", Dist {"1/6 * T * (T^6 - 1) / (T - 1)"}, "n");
    cout << "l2:\t" << l2 << endl;
    l3 = l2.update("m", "m+t");
    cout << "l3:\t" << l3 << endl;
    l4 = l3.update("n", "0");
    cout << "l4:\t" << l4 << endl;
    l5 = l4.update("t", "0");
    cout << "l5:\t" << l5 << endl;
    Dist res1 = l5 + l0.filterLeq("n", "0");
    cout << "res1:\t" << res1 << endl;

    /*
     * guard + body + invariant
     *
     0:
     1: if (n > 0) {
     2:   tmp := unif(1,6);
     3:   m := m + tmp;
     4:   tmp := 0;
     5:   n := n-1;
     6:   if (0 < n){
     7:     tmp := iid(unif(1,6), n)
     8:     m := m+tmp
     9:     n := 0
    10:     tmp := 0
    11:   }
    12:}
     *
     */
    cout << "Modified Inv" << endl;
    l0 = Dist {"1/(1 - p*n)", "p"} *  Dist {"1/(1 - q*m)", "q"} * Dist {"1/(1 - r*t)", "r"};
    l1 = l0.filterGreater("n", "0");
    l2 = l1.update("t", "0") * Dist {"1/6 * T * (T^6 - 1) / (T - 1)"};
    l3 = l2.update("m", "m+t");
    l4 = l3.update("t", "0");
    l5 = l4.update("n", "n-1");
    l6 = l5.filterGreater("n", "0");
    l7 = l6.updateIid("t", Dist {"1/6 * T * (T^6 - 1) / (T - 1)"}, "n");
    l8 = l7.update("m", "m+t");
    l9 = l8.update("n", "0");
    l10 = l9.update("t", "0");
    l11 = l10 + l5.filterLeq("n", "0");
    l12 = l11 + l0.filterLeq("n", "0");
    Dist res2 = l12;

    ASSERT_EQ(res1, res2);
}

TEST(EquivalenceCheck, dep_bern){
    Dist::resetSymbolCache();
    Dist l0, l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12;
    
    /*
     0:
     1: while( c > 0 ){
     2:    {m := m+1} [1/2] {n:= n+1}
     3:    c := c-1
     4:    tmp := 0
     5: }
     */
    
    /*
     0: Invariant
     1: if (c > 0){
     2:    tmp := iid(bernoulli(1/2), c)
     3:    m := m + tmp
     4:    n := n + (c-tmp)
     5:    c := 0
     6:    tmp := 0
     7: } else {skip}
     */
    
    cout << "Invariant" << endl;
    l0 = Dist {"1/(1 - p*n)", "p"} *  Dist {"1/(1 - q*m)", "q"} * Dist {"1/(1 - r*c)", "r"} * Dist {"1/(1-s*t)", "s"};
    l1 = l0.filterGreater("c", "0");
    l2 = l1.updateIid("t", Dist {"1/2 * T + 1/2"}, "c");
    l3 = l2.update("m", "m+t");
    l4 = l3.update("n", "n+c-t");
    l5 = l4.update("c", "0");
    l6 = l5.update("t", "0");
    l7 = l6 + (l0 - l1);
    Dist res1 = l7;
    
    cout << res1 << endl;
    
    /*
     0: Guard + Body + Invariant
     1: if (c > 0){
     2:    {m := m+1} [1/2] {n:= n+1};
     3:    c := c-1;
     4:    tmp := 0;
     5:    if (c > 0){
     6:       tmp := iid(bernoulli(1/2), c)
     7:       m := m + tmp
     8:       n := n + (c-tmp)
     9:       c := 0
     10:      tmp := 0
     11:    } else {skip}
     12: }
     */
    
    cout << "Modified Invariant" << endl;
    l0 = Dist {"1/(1 - p*n)", "p"} *  Dist {"1/(1 - q*m)", "q"} * Dist {"1/(1 - r*c)", "r"} * Dist {"1/(1-s*t)", "s"};
    l1 = l0.filterGreater("c", "0");
    l2 = "1/2" * l1.update("m", "m+1") + "1/2" * l1.update("n", "n+1");
    l3 = l2.update("c", "c-1");
    l4 = l3.update("t", "0");
    l5 = l4.filterGreater("c", "0");
    l6 = l5.updateIid("t", Dist {"1/2 * T + 1/2"}, "c");
    l7 = l6.update("m", "m+t");
    l8 = l7.update("n", "n+c-t");
    l9 = l8.update("c", "0");
    l10 = l9.update("t", "0");
    l11 = l10 + l4.filterLeq("c", "0");
    l12 = l11 + l0.filterLeq("c", "0");
    Dist res2 = l12;
    
    cout << res2 << endl;
    
    ASSERT_EQ(res1, res2);
}
