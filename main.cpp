// компілятор - MSVC
// варіант 1

#include <chrono>
#include <execution>
#include <format>
#include <functional>
#include <iostream>
#include <latch>
#include <map>
#include <random>
#include <vector>

bool MyPredicate(const int number) {
    return (2 * atan(tanh(number / 2.0))) > 0;
}

template <typename Function, typename... Arguments>
std::chrono::duration<double, std::milli> MeasureFunctionExecutionTime(
    Function &&function, const Arguments &...arguments) {
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

    std::uniform_int_distribution<std::mt19937::result_type> distribution(2, 8);

    for (int i = 0; i < sequence_size; i++) {
        random_sequence.push_back(distribution(mersenne_twister_engine));
    }

    return random_sequence;
}

template <typename VectorElementType>
std::vector<std::vector<int>> SplitSequenceIntoChunks(
    const std::vector<VectorElementType> &sequence, int chunk_count) {
    std::vector<std::vector<int>> chunks(chunk_count);
    const int chunk_size = sequence.size() / chunk_count;
    const int remainder = sequence.size() % chunk_count;

    int start = 0;
    for (int i = 0; i < chunk_count; ++i) {
        const int end = start + chunk_size + (i < remainder ? 1 : 0);
        chunks[i] =
            std::vector(sequence.begin() + start, sequence.begin() + end);
        start = end;
    }

    return chunks;
}

void TestAllOfUsingMultipleThreads(
    const std::vector<std::vector<int>> &sequence_chunks) {
    std::vector<std::jthread> threads;
    std::latch work_done(sequence_chunks.size());

    for (const auto &sequence_chunk : sequence_chunks) {
        threads.emplace_back([&sequence_chunk, &work_done]() {
            auto _ = std::ranges::all_of(sequence_chunk, [](const int number) {
                return MyPredicate(number);
            });
            work_done.count_down();
        });
    }

    work_done.wait();
}

void Part1(const std::vector<std::vector<int>> &sequences) {
    std::cout << "<start> No policy\n";

    for (const auto &sequence : sequences) {
        std::cout << std::format("sequence size = {} ; ", sequence.size())
                  << "time = " << MeasureFunctionExecutionTime([&sequence]() {
                         auto _ = std::ranges::all_of(sequence, MyPredicate);
                     })
                  << std::endl;
    }
    std::cout << "<end>\n\n" << std::endl;
}

void Part2(const std::vector<std::vector<int>> &sequences) {
    std::cout << "<start> Sequenced policy\n";
    for (const auto &sequence : sequences) {
        std::cout << std::format("sequence size = {} ; ", sequence.size())
                  << "time = " << MeasureFunctionExecutionTime([&sequence]() {
                         auto _ =
                             std::all_of(std::execution::seq, sequence.begin(),
                                         sequence.end(), MyPredicate);
                     })
                  << std::endl;
    }
    std::cout << "<end>\n\n" << std::endl;

    std::cout << "<start> Parallel policy\n";
    for (const auto &sequence : sequences) {
        std::cout << std::format("sequence size = {} ; ", sequence.size())
                  << "time = " << MeasureFunctionExecutionTime([&sequence]() {
                         auto _ =
                             std::all_of(std::execution::par, sequence.begin(),
                                         sequence.end(), MyPredicate);
                     })
                  << std::endl;
    }
    std::cout << "<end>\n\n" << std::endl;

    std::cout << "<start> Parallel unsequenced policy\n";
    for (const auto &sequence : sequences) {
        std::cout << std::format("sequence size = {} ; ", sequence.size())
                  << "time = " << MeasureFunctionExecutionTime([&sequence]() {
                         auto _ = std::all_of(std::execution::par_unseq,
                                              sequence.begin(), sequence.end(),
                                              MyPredicate);
                     })
                  << std::endl;
    }
    std::cout << "<end>\n\n" << std::endl;

    std::cout << "<start> Unsequenced policy\n";
    for (const auto &sequence : sequences) {
        std::cout << std::format("sequence size = {} ; ", sequence.size())
                  << "time = " << MeasureFunctionExecutionTime([&sequence]() {
                         auto _ = std::all_of(std::execution::unseq,
                                              sequence.begin(), sequence.end(),
                                              MyPredicate);
                     })
                  << std::endl;
    }
    std::cout << "<end>\n\n" << std::endl;
}

struct ChunkedSequence {
    std::vector<std::vector<int>> chunks;
    std::size_t sequence_size;
};

void Part3(const std::vector<std::vector<int>> &sequences) {
    std::map<int, std::vector<ChunkedSequence>>
        chunk_count_to_chunked_sequences_map;

    for (const auto &sequence : sequences) {
        for (int k = 2; k < 50; k += 2) {
            if (!chunk_count_to_chunked_sequences_map.contains(k)) {
                chunk_count_to_chunked_sequences_map[k] =
                    std::vector<ChunkedSequence>();
            }

            ChunkedSequence split_sequence;
            split_sequence.chunks = SplitSequenceIntoChunks(sequence, k);
            split_sequence.sequence_size = sequence.size();

            chunk_count_to_chunked_sequences_map[k].push_back(split_sequence);
        }
    }

    for (auto const &[chunk_count, chunked_sequences] :
         chunk_count_to_chunked_sequences_map) {
        std::cout << std::format("<start> Multithreaded k = {} \n",
                                 chunk_count);

        for (const auto &chunked_sequence : chunked_sequences) {
            std::cout << std::format("sequence size = {} ; ",
                                     chunked_sequence.sequence_size)
                      << "time = "
                      << MeasureFunctionExecutionTime(
                             TestAllOfUsingMultipleThreads,
                             chunked_sequence.chunks)
                      << std::endl;
        }

        std::cout << "<end>\n\n" << std::endl;
    }
}

int main() {
    std::random_device random_device;
    std::mt19937 mersenne_twister_engine(random_device());

    constexpr int sequence_count = 5;
    std::vector<std::vector<int>> sequences;
    int starting_sequence_size = 100;
    sequences.reserve(sequence_count);

    for (int i = 0; i < sequence_count; i++) {
        sequences.push_back(GenerateRandomSequence(starting_sequence_size,
                                                   mersenne_twister_engine));
        starting_sequence_size *= 10;
    }

    Part1(sequences);
    Part2(sequences);
    Part3(sequences);

    return 0;
}