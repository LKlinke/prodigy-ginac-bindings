//
// Created by Tobias Winkler on 01.11.21.
//

#include "gtest/gtest.h"
#include "Dist.h"
#include "DistributionFactory.h"

using namespace std;

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
    // l0 = DistributionFactory::geometric("x", "p") * DistributionFactory::geometric("c", "q");
    l0 = Dist {"1/(1 - p*x)", "p"} * Dist {"1/ (1 - q*c)", "q"};
    l1 = l0.updateIidSampling("tmp", DistributionFactory::geometric("TMP", "1/2"), "x");
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
    l0 = Dist {"1/(1 - p*x)", "p"} * Dist {"1/ (1 - q*c)", "q"};
    l1 = l0.filterGreater("x", "0");
    l2 = "1/2" * l1.update("x", "x - 1") + "1/2" * l1.update("c", "c + 1");
    l3 = l2.updateIidSampling("tmp", DistributionFactory::geometric("TMP", "1/2"), "x");
    l4 = l3.update("c", "c + tmp");
    l5 = l4.update("tmp", "0");
    l6 = l5.update("x", "0");
    l7 = l6 + (l0 - l1);
    Dist res2 = l7;

    ASSERT_EQ(res1, res2);
}