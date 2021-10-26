import prodigy

if __name__ == "__main__":
    dist = prodigy.geometric("x", "1/2")
    print(dist.E("x"))
    print(dist.hasParameter("q"))
    dist = dist * "q"
    print(dist.hasParameter("q"))
    print(dist.E("x"))

    dist = prodigy.geometric("x", "1/2")
    dist = dist.update("y", "x + 5")
    print(dist.E("x"))
