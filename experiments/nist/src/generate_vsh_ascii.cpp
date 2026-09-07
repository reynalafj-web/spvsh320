#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <openssl/bn.h>
#include <iomanip>

// Fungsi untuk mengecek bilangan prima (untuk membangkitkan 256 prima pertama)
bool is_prime(int n) {
    if (n < 2) return false;
    for (int i = 2; i * i <= n; i++) {
        if (n % i == 0) return false;
    }
    return true;
}

int main() {
    const uint64_t TARGET_BITS = 1000000000; // 1 Miliar Bit
    const int DIGEST_BITS = 256;
    const uint64_t TARGET_ITERATIONS = TARGET_BITS / DIGEST_BITS; // 3.906.250 iterasi
    const std::string FILENAME = "vsh_final_1Gb.txt";

    std::cout << "[INFO] Membangkitkan 256 bilangan prima pertama..." << std::endl;
    std::vector<int> primes;
    int num = 2;
    while (primes.size() < 256) {
        if (is_prime(num)) primes.push_back(num);
        num++;
    }

    // Inisialisasi OpenSSL BIGNUM
    BN_CTX *ctx = BN_CTX_new();
    BIGNUM *n = BN_new();
    BIGNUM *x = BN_new();
    BIGNUM *prod = BN_new();
    BIGNUM *tmp = BN_new();

    // Modulus 256-bit (Gunakan nilai yang sama dengan benchmark sebelumnya)
    BN_hex2bn(&n, "C4B36F8B6774640822BEF3DEB3E2F49BA5EDAC11812831E63C1B609E839B6B4D");

    std::ofstream outfile(FILENAME);
    if (!outfile.is_open()) {
        std::cerr << "Gagal membuat file!" << std::endl;
        return 1;
    }

    std::cout << "[INFO] Memulai generasi 1 Miliar bit VSH-256 (Counter-mode)..." << std::endl;

    for (uint64_t counter = 0; counter < TARGET_ITERATIONS; counter++) {
        // VSH Core: x = 1
        BN_one(x);

        // Payload 8-byte (counter) diproses sebagai 64 bit pertama
        // Di VSH asli, bit ke-j menentukan apakah prima ke-j dikalikan
        for (int j = 0; j < 64; j++) {
            if ((counter >> j) & 1) {
                BN_set_word(tmp, primes[j]);
                BN_mul(x, x, tmp, ctx);
            }
        }
        
        // Finalisasi modulo
        BN_mod(x, x, n, ctx);

        // Konversi hasil BIGNUM ke bit string ASCII
        unsigned char bin_out[32];
        BN_bn2binpad(x, bin_out, 32);
        
        std::string bit_string = "";
        for (int b = 0; b < 32; b++) {
            for (int i = 7; i >= 0; i--) {
                bit_string += ((bin_out[b] >> i) & 1) ? '1' : '0';
            }
        }

        outfile << bit_string;

        if (counter % 100000 == 0) {
            std::cout << "\rProgress: " << (counter * 256) / 1000000 << " / 1000 Mb" << std::flush;
        }
    }

    outfile.close();
    std::cout << "\n[SUKSES] File " << FILENAME << " berhasil dibuat." << std::endl;

    // Free memory
    BN_free(n); BN_free(x); BN_free(prod); BN_free(tmp);
    BN_CTX_free(ctx);

    return 0;
}
