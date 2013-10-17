from Ivy.benchmark import Benchmark

if __name__ == '__main__':
    ans = [291, 393, 93, 33, 1, 3, 4,6,0]
    pred = [9,32,1,99,93, 939, 10]
    bench = Benchmark(answer=ans, predict=pred)
    
    print bench.recall()
    print bench.precision()
    print bench.f_measure()
