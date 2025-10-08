// компілятор - MSVC
// варіант 1

#include <chrono>
#include <execution>
#include <functional>
#include <iostream>
#include <latch>
#include <random>

bool IsEven(const int number) { return number % 2 == 0; }

void TestAllOfWithoutPolicy(const std::vector<int> &sequence) {
    auto _ = std::all_of(sequence.begin(), sequence.end(),
                         [](const int number) { return IsEven(number); });
}

void TestAllOfWithSequencedPolicy(const std::vector<int> &sequence) {
    auto _ = std::all_of(std::execution::seq, sequence.begin(), sequence.end(),
                         [](const int number) { return IsEven(number); });
}

void TestAllOfWithParallelPolicy(const std::vector<int> &sequence) {
    auto _ = std::all_of(std::execution::par, sequence.begin(), sequence.end(),
                         [](const int number) { return IsEven(number); });
}

void TestAllOfWithParallelUnsequencedPolicy(const std::vector<int> &sequence) {
    auto _ =
        std::all_of(std::execution::par_unseq, sequence.begin(), sequence.end(),
                    [](const int number) { return IsEven(number); });
}

void TestAllOfWithUnsequencedPolicy(const std::vector<int> &sequence) {
    auto _ =
        std::all_of(std::execution::unseq, sequence.begin(), sequence.end(),
                    [](const int number) { return IsEven(number); });
}

template <typename Function, typename... Arguments>
std::chrono::duration<double, std::milli> MeasureFunctionExecutionTime(
    const Function &&function, const Arguments &...arguments) {
    const auto time_point_before_function_call =
        std::chrono::high_resolution_clock::now();
    function(arguments...);
    const auto time_point_after_function_call =
        std::chrono::high_resolution_clock::now();
    return time_point_after_function_call - time_point_before_function_call;
}

std::vector<int> GenerateRandomSequence(const std::size_t sequence_size,
                                        std::mt19937 &mersenne_twister_engine) {
    std::vector<int> random_sequence;
    random_sequence.reserve(sequence_size);

    std::uniform_int_distribution<std::mt19937::result_type> distribution(1, 6);

    for (int i = 0; i < sequence_size; i++) {
        random_sequence.push_back(distribution(mersenne_twister_engine));
    }

    return random_sequence;
}

template <typename VectorElementType>
std::vector<std::vector<int>> SplitVectorIntoKVectors(
    const std::vector<VectorElementType> &vector_to_split, int k) {
    std::vector<std::vector<int>> chunks(k);
    const int chunk_size = vector_to_split.size() / k;
    const int remainder = vector_to_split.size() % k;

    int start = 0;
    for (int i = 0; i < k; ++i) {
        const int end = start + chunk_size + (i < remainder ? 1 : 0);
        chunks[i] = std::vector(vector_to_split.begin() + start,
                                vector_to_split.begin() + end);
        start = end;
    }

    return chunks;
}

void TestAllOfUsingMultipleThreads(
    const std::vector<std::vector<int>> &sequence_chunks) {
    std::vector<std::jthread> threads;
    std::latch work_done(sequence_chunks.size());
    for (int i = 0; i < sequence_chunks.size(); i++) {
        threads.emplace_back([sequence_chunks, i, &work_done]() {
            auto _ = std::ranges::all_of(
                sequence_chunks[i],
                [](const int number) { return IsEven(number); });
            work_done.count_down();
        });
    }

    work_done.wait();
}

int main() {
    std::random_device random_device;
    std::mt19937 mersenne_twister_engine(random_device());

    const auto seq1 = GenerateRandomSequence(10000, mersenne_twister_engine);
    const auto seq2 = GenerateRandomSequence(100000, mersenne_twister_engine);
    const auto seq3 = GenerateRandomSequence(1000000, mersenne_twister_engine);

    // part 1
    std::cout << "<start> No policy\n";
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithoutPolicy, seq1)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithoutPolicy, seq2)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithoutPolicy, seq3)
              << std::endl;
    std::cout << "<end>\n\n" << std::endl;

    // part 2
    std::cout << "<start> Sequenced policy\n";
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithSequencedPolicy,
                                              seq1)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithSequencedPolicy,
                                              seq2)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithSequencedPolicy,
                                              seq3)
              << std::endl;
    std::cout << "<end>\n\n" << std::endl;

    std::cout << "<start> Parallel policy\n";
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithParallelPolicy, seq1)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithParallelPolicy, seq2)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithParallelPolicy, seq3)
              << std::endl;
    std::cout << "<end>\n\n" << std::endl;

    std::cout << "<start> Parallel unsequenced policy\n";
    std::cout << MeasureFunctionExecutionTime(
                     TestAllOfWithParallelUnsequencedPolicy, seq1)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(
                     TestAllOfWithParallelUnsequencedPolicy, seq2)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(
                     TestAllOfWithParallelUnsequencedPolicy, seq3)
              << std::endl;
    std::cout << "<end>\n\n" << std::endl;

    std::cout << "<start> Unsequenced policy\n";
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithUnsequencedPolicy,
                                              seq1)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithUnsequencedPolicy,
                                              seq2)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfWithUnsequencedPolicy,
                                              seq3)
              << std::endl;
    std::cout << "<end>\n\n" << std::endl;

    // part 3
    const int k = 100;
    const auto seq1_split = SplitVectorIntoKVectors(seq1, k);
    const auto seq2_split = SplitVectorIntoKVectors(seq2, k);
    const auto seq3_split = SplitVectorIntoKVectors(seq3, k);

    std::cout << "<start> Multithreaded \n";
    std::cout << MeasureFunctionExecutionTime(TestAllOfUsingMultipleThreads,
                                              seq1_split)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfUsingMultipleThreads,
                                              seq2_split)
              << std::endl;
    std::cout << MeasureFunctionExecutionTime(TestAllOfUsingMultipleThreads,
                                              seq3_split)
              << std::endl;
    std::cout << "<end>\n\n" << std::endl;

    return 0;
}