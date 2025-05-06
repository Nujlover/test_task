#include <iostream>
#include <vector>
#include <complex>
#include <random>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <limits>

using namespace std;

class QAMModulator {
    int modulation_order;
    vector<complex<double>> constellation;

public:
    QAMModulator(int mod_order) : modulation_order(mod_order) {
        int bits_per_symbol = static_cast<int>(log2(mod_order)); 
        int side = static_cast<int>(sqrt(mod_order));            

        for (int i = 0; i < mod_order; ++i) {
            int x = (i % side) - side / 2;
            int y = (i / side) - side / 2;
            constellation.emplace_back(2 * x + 1, 2 * y + 1);
        }

        double norm = sqrt((mod_order - 1)*2/ 3.0);
        for (auto& c : constellation) c /= norm;

    }

    vector<complex<double>> modulate(const vector<bool>& bits) {
        vector<complex<double>> symbols;
        int bits_per_symbol = static_cast<int>(log2(modulation_order));

        if (bits.size() % bits_per_symbol != 0) {
            throw invalid_argument("Количество бит не кратно " + to_string(bits_per_symbol));
        }

        for (size_t i = 0; i < bits.size(); i += bits_per_symbol) {
            int index = 0;
            for (int j = 0; j < bits_per_symbol; ++j) {
                index = (index << 1) | static_cast<int>(bits[i + j]); 
            }
            symbols.push_back(constellation[index]);
        }
        return symbols;
    }
};


class NoiseAdder {
    normal_distribution<double> dist;
    mt19937 generator;

public:
    NoiseAdder(double variance) : dist(0.0, sqrt(variance / 2)), generator(random_device{}()) {}

    vector<complex<double>> addNoise(const vector<complex<double>>& symbols) {
        vector<complex<double>> noisy(symbols);
        for (auto& s : noisy) {
            s += complex<double>(dist(generator), dist(generator));
        }
        return noisy;
    }
};


class QAMDemodulator {
    int modulation_order;
    vector<complex<double>> constellation;

public:
    QAMDemodulator(int mod_order) : modulation_order(mod_order) {
        int side = static_cast<int>(sqrt(mod_order));
        for (int i = 0; i < mod_order; ++i) {
            int x = (i % side) - side / 2;
            int y = (i / side) - side / 2;
            constellation.emplace_back(2 * x + 1, 2 * y + 1);
        }
        double norm = sqrt((mod_order - 1)*2/ 3.0);
        for (auto& c : constellation) c /= norm;
    }

    vector<bool> demodulate(const vector<complex<double>>& symbols) {
        vector<bool> bits;
        int bits_per_symbol = static_cast<int>(log2(modulation_order));

        for (const auto& s : symbols) {
            double min_dist = numeric_limits<double>::max();
            int best_index = 0;

            for (int i = 0; i < modulation_order; ++i) {
                double d = norm(s - constellation[i]);
                if (d < min_dist) {
                    min_dist = d;
                    best_index = i;
                }
            }

            for (int j = bits_per_symbol - 1; j >= 0; --j) {
                bits.push_back((best_index >> j) & 1);
            }
        }
        return bits;
    }
};


vector<bool> generate_bits(size_t num_bits) {
    random_device rd;
    mt19937 gen(rd());
    bernoulli_distribution d(0.5);

    vector<bool> bits(num_bits);
    for (size_t i = 0; i < num_bits; ++i) { 
        bits[i] = d(gen);
    }
    return bits;
}


double calculate_ber(const vector<bool>& tx, const vector<bool>& rx) {
    if (tx.size() != rx.size()) throw invalid_argument("Size mismatch");

    size_t errors = 0;
    for (size_t i = 0; i < tx.size(); ++i) {
        errors += (tx[i] != rx[i]);
    }
    return static_cast<double>(errors) / tx.size();
}

int main() {
    const int MODULATION = 16; //4 - QPSK; 16 - 16-QAM; 64 - 64_QAM
    const vector<double> variances = { 0.00316, 0.01,0.0316, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0 };

    const size_t NUM_BITS = log2(MODULATION) * 1000000;

    auto tx_bits = generate_bits(NUM_BITS);
    QAMModulator modulator(MODULATION);
    QAMDemodulator demodulator(MODULATION);

    vector<double> bers;
    vector<double> snrs;  

    for (auto var : variances) {
        auto symbols = modulator.modulate(tx_bits);
        NoiseAdder noise_adder(var);
        auto noisy_symbols = noise_adder.addNoise(symbols);
        auto rx_bits = demodulator.demodulate(noisy_symbols);
        double ber = calculate_ber(tx_bits, rx_bits);
        bers.push_back(ber);

        double snr_db = -10.0 * log10(var); 
        snrs.push_back(snr_db);

        cout << "dispersion: " << var
            << " | SNR (dB): " << snr_db
            << " | BER: " << ber << endl;
    }

    ofstream out("ber_results.txt");
    for (size_t i = 0; i < variances.size(); ++i) {
        out << variances[i] << " " << snrs[i] << " " << bers[i] << endl;
    }

    return 0;
}