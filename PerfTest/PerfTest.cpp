// PerfTest.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <chrono>
#include <thread>
#include <windows.h>
struct xoshiro256p_state
{
    uint64_t s[4];
};

uint64_t rol64(uint64_t x, int k)
{
    return (x << k) | (x >> (64 - k));
}

uint64_t xoshiro256p(struct xoshiro256p_state* state)
{
    uint64_t* s = state->s;
    uint64_t const result = s[0] + s[3];
    uint64_t const t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = rol64(s[3], 45);

    return result;
}

uint64_t lcg(uint64_t state) {
    return state * 2862933555777941757 + 3037000493;
}

uint64_t splitmix(uint64_t state) {
    uint64_t z = (state += 0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
    z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
    return z ^ (z >> 31);
}

xoshiro256p_state xoshiro256p_init(uint64_t seed)
{
	struct xoshiro256p_state state;
	state.s[0] = splitmix(seed);
	state.s[1] = splitmix(state.s[0]);
	state.s[2] = splitmix(state.s[1]);
	state.s[3] = splitmix(state.s[2]);

    for (int i = 0; i < 100; i++)
    {
		xoshiro256p(&state);
	}

	return state;
}

// perform mem copy test, return speed in gigabytes/second
uint64_t performCopyTest(size_t count, int threadId, size_t* startval, double* scores)
{
    size_t* array1 = new size_t[count];

    // initialize array1

    struct xoshiro256p_state state = xoshiro256p_init(time(NULL) + threadId);

    for (size_t i = 0; i < count; i++)
    {
        array1[i] = xoshiro256p(&state);
    }

    size_t iters = 2.0 / ((double)count / 1024.0 / 1024.0 / 1024.0);

    while (startval[0] != 0x1234)
    {
		// wait until start is 0x1234
	}

    auto start = std::chrono::high_resolution_clock::now();

    size_t* src = array1;

    uint64_t sum = 0;
    for (size_t i = 0; i < iters; i++)
    {
        //std::copy(src, src + count, dst);

        for (size_t j = 0; j < count/8; j++)
        {
            const size_t a0 = src[j * 8 + 0];
            const size_t a1 = src[j * 8 + 1];
            const size_t a2 = src[j * 8 + 2];
            const size_t a3 = src[j * 8 + 3];
            const size_t a4 = src[j * 8 + 4];
            const size_t a5 = src[j * 8 + 5];
            const size_t a6 = src[j * 8 + 6];
            const size_t a7 = src[j * 8 + 7];

            sum += a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7;
        }
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    double seconds = elapsed.count() / iters;
    double speed = count * sizeof(size_t) / seconds / 1024 / 1024 / 1024;

    scores[threadId] = speed;

    delete[] array1;

    return sum;
}


uint64_t xoshirogen(size_t count, int threadId, size_t* startval, double* scores)
{

    struct xoshiro256p_state state = xoshiro256p_init(time(NULL) + threadId);

    size_t iters = 100000000;

    while (startval[0] != 0x1234)
    {
        // wait until start is 0x1234
    }

    auto start = std::chrono::high_resolution_clock::now();

    uint64_t sum = 0;
    for (size_t i = 0; i < iters; i++)
    {
		sum += xoshiro256p(&state);
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    double seconds = elapsed.count();
	size_t generated = iters * sizeof(uint64_t);

	std::cout << sum << "\n";

	std::cout << "thread " << threadId << " generated " << generated << " bytes in " << seconds << " seconds\n";

	scores[threadId] = generated / seconds / 1024 / 1024 / 1024;

    return sum;
}

uint64_t memRandAccess(size_t count, int threadId, size_t* startval, double* scores)
{
    size_t* array1 = new size_t[count];

    // initialize array1

    struct xoshiro256p_state state = xoshiro256p_init(time(NULL) + threadId);

    for (size_t i = 0; i < count; i++)
    {
		array1[i] = xoshiro256p(&state) % count;
    }

    size_t iters = count * 1000;

    while (startval[0] != 0x1234)
    {
        // wait until start is 0x1234
    }

    size_t next = 0;
    for (size_t i = 0; i < iters; i++)
    {
        next = array1[next];
    }

    auto start = std::chrono::high_resolution_clock::now();

    uint64_t sum = 0;
    for (size_t i = 0; i < iters; i++)
    {
        next = array1[next];

		sum += next;
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    double seconds = elapsed.count() / iters;

	scores[threadId] = seconds * 1000000000;

    delete[] array1;

    return sum;
}

std::atomic<size_t> testatomic(0);

void c2ctest(int threadId, int otherThreadId, size_t* startval, double* scores)
{
    while (startval[0] != 0x1234)
    {
        // wait until start is 0x1234
    }

    size_t iters = 1500000;

    auto start = std::chrono::high_resolution_clock::now();

    for (size_t i = 0; i < iters; i++)
    {
        size_t expected = otherThreadId;
        while (!testatomic.compare_exchange_weak(expected, threadId))
        {
			expected = otherThreadId;
		}
    }

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    double seconds = elapsed.count();

    double speed = (iters * 2) / seconds;

    if (scores == nullptr)
    {
        return;
    }

    int hwThreads = std::thread::hardware_concurrency();

    scores[threadId * hwThreads + otherThreadId] = 1 / speed * 1000000000;
}

void performc2c()
{
    int hwThreads = std::thread::hardware_concurrency();

    std::cout << "Hardware threads: " << hwThreads << std::endl;

    size_t* testval = new size_t[1];
    double *scores = new double[hwThreads * hwThreads];

    for (size_t i = 0; i < hwThreads; i++)
    {
        for (size_t j = 0; j < hwThreads; j++)
        {
            if (i == j)
            {
                scores[i * hwThreads + j] = 0;
            }

            testval[0] = j;
            size_t* start = new size_t[1];

            testatomic.store(j);

            std::thread threadA(c2ctest, i, j, start, scores);
            std::thread threadB(c2ctest, j, i, start, nullptr);

            SetThreadAffinityMask(threadA.native_handle(), 1 << i);
            SetThreadAffinityMask(threadB.native_handle(), 1 << j);

            // theads are spinning, start the test
            start[0] = 0x1234;

            threadA.join();
            threadB.join();

            std::cout << "thread " << i << " to " << j << " latency: " << scores[i * hwThreads + j] << " ns\n";

            delete[] start;
        }
	}

    for (size_t i = 0; i < hwThreads; i++)
    {
        for (size_t j = 0; j < hwThreads; j++)
        {
			std::cout << scores[i * hwThreads + j] << ", ";
		}
		std::cout << ";\n";
    }

    delete[] scores;
}

int main()
{
    std::cout << "Hello World!\n";

    size_t threads = 32;

    //performc2c();

    //return 0;
    for (size_t s = 1; s < 24; s++)
    {
        size_t count = pow(2,s) * 128 / sizeof(size_t);

        std::thread* thread = new std::thread[threads];

        double* scores = new double[threads];

        size_t* start = new size_t[1];

        for (size_t i = 0; i < threads; i++)
        {
            thread[i] = std::thread(xoshirogen, count, i,start, scores);
            SetThreadAffinityMask(thread[i].native_handle(), 1 << i);
        }

        start[0] = 0x1234;

        for (size_t i = 0; i < threads; i++)
        {
            thread[i].join();
        }

        double speed = 0;

        for (size_t i = 0; i < threads; i++)
        {
            speed += scores[i];
        }

        std::cout << "work size: " << count * sizeof(size_t) /1024 << " KB, speed: " << speed << " ns\n";
    }
}