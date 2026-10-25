// Exact comparisons of parallel signed floors and guard certificates.
#define main weighted_parallel_program_main
#include "../decimal_extract_weighted_machin_parallel.cpp"
#undef main

int main() {
    try {
        u64 checks=0,certificates=0;
        for(u64 N:{1ULL,2ULL,17ULL,100ULL,257ULL,1024ULL,4097ULL}) {
            for(u64 guard:{1ULL,17ULL,63ULL,64ULL,127ULL,192ULL,257ULL}) {
                machin::Stats serial; serial.block_length=31;
                const auto expected=machin::serial_attempt(N,guard,serial);
                for(u64 threads:{2ULL,3ULL,8ULL}) {
                    machin::Stats actual; actual.block_length=31; actual.threads=threads;
                    const auto residue=machin::parallel_attempt(N,guard,actual);
                    if(residue!=expected || actual.terms!=serial.terms ||
                       actual.word_errors!=serial.word_errors || actual.core_errors!=serial.core_errors ||
                       actual.core_terms!=serial.core_terms || actual.weighted.blocks!=serial.weighted.blocks ||
                       actual.weighted.real_terms!=serial.weighted.real_terms)
                        throw std::runtime_error("parallel attempt differs from exact serial floors or counts");
                    ++checks;
                }
            }
        }
        for(u64 N:{1ULL,17ULL,257ULL,4097ULL}) {
            machin::Stats serial,parallel;
            serial.threads=1; parallel.threads=8;
            const auto expected=machin::extract(N,1,false,serial);
            const auto actual=machin::extract(N,1,false,parallel);
            if(actual!=expected || serial.guard!=parallel.guard || serial.passes!=parallel.passes ||
               serial.word_errors!=parallel.word_errors || serial.core_errors!=parallel.core_errors)
                throw std::runtime_error("parallel guard retries or certificate differ");
            if(serial.passes<=1) throw std::runtime_error("retry test failed to exercise adaptive guards");
            ++certificates;
        }
        lcm_invariant_checks=true;
        machin::Stats fallback; fallback.threads=8;
        const auto a=machin::attempt(257,192,fallback);
        if(fallback.worker_jobs) throw std::runtime_error("invariant diagnostics were parallelized");
        machin::Stats serial;
        if(a!=machin::serial_attempt(257,192,serial)) throw std::runtime_error("diagnostic fallback mismatch");
        std::cout<<"{\"parallel_checks\":\"PASS\",\"exact_attempts\":"<<checks
                 <<",\"adaptive_certificates\":"<<certificates
                 <<",\"hardware_threads\":"<<std::thread::hardware_concurrency()<<"}\n";
        return 0;
    } catch(const std::exception& error) {std::cerr<<error.what()<<'\n'; return 1;}
}
