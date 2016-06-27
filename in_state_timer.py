import timeit

if __name__ == "__main__":

    setup = """
from RegionTester.region_tester import InUSState
California = InUSState('')
California.get_state_shape('California')
"""
    t = timeit.timeit('California.point_inside(-120, 38)',
                      setup=setup,
                      number=10)
    print t
    t = timeit.timeit('California.point_inside(-20, 45)',
                      setup=setup,
                      number=10)
    print "time: ", t
