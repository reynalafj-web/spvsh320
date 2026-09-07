# Cara mengunggah paket ini ke GitHub dan Zenodo

PeerJ Computer Science tidak menerima tautan GitHub yang bisa berubah sebagai satu-satunya artefak. Alur yang diterima reviewer: **GitHub (kerja) + tag beku + arsip Zenodo ber-DOI**.

Ganti tiga placeholder sebelum submit:

- `YOUR_GITHUB_USER`
- `spvsh320` (nama repo; boleh diganti)
- email penulis korespondensi

---

## A. Siapkan akun

1. GitHub: https://github.com/signup
2. Zenodo: https://zenodo.org/signup
3. Login Zenodo memakai GitHub: https://zenodo.org/account/settings/github/
4. Di halaman itu, nyalakan sakelar untuk repositori `spvsh320` setelah repo ada.

---

## B. Buat repositori GitHub

Di komputer lokal (bukan wajib di server eksperimen):

```bash
# unpack zip paket ini, lalu
cd spvsh320-repro

git init
git add .
git commit -m "Add paper-v1 reproducibility snapshot"

git branch -M main
git remote add origin https://github.com/YOUR_GITHUB_USER/spvsh320.git
git push -u origin main
```

Di situs GitHub:

1. Create repository → nama `spvsh320` → Public.
2. Jangan centang “Add README” (paket ini sudah punya README).
3. Settings → Pages tidak diperlukan.
4. About → Add topics: `cryptography`, `hash-function`, `reproducibility`, `peerj`.
5. Releases → jangan dibuat dulu sampai langkah C selesai.

Visibility: **Public** untuk publikasi. Untuk review saja, repo Private + tautan reviewer juga diterima sementara, tetapi DOI Zenodo publik tetap harus ada sebelum accepted paper terbit.

---

## C. Bekukan versi yang sama dengan PDF

```bash
git tag -a paper-v1 -m "Snapshot matching the submitted PeerJ manuscript"
git push origin paper-v1
```

Lalu di GitHub: Releases → Draft a new release

- Tag: `paper-v1`
- Title: `SpVSH-320 reproducibility package paper-v1`
- Description: satu paragraf dari README (target, kandidat, rasio utama)
- Attach binaries? Tidak perlu. Source + CSV sudah di tag.

Jangan menambah file ke `main` setelah tag jika file itu dikutip di PDF. Perubahan sesudah submit = tag baru (`paper-v1.1`) dan versi Zenodo baru.

---

## D. Arsip ke Zenodo

Setelah sakelar GitHub–Zenodo aktif:

1. Publish release `paper-v1` di GitHub.
2. Zenodo akan membuat draft record otomatis (kadang 2–10 menit).
3. Buka https://zenodo.org/account/settings/github/ lalu klik record tersebut.
4. Lengkapi metadata sebelum Publish:

| Field | Isi yang disarankan |
|---|---|
| Title | SpVSH-320 reproducibility package |
| Creators | Reynaldi Alfajri; Farah Afianti |
| Affiliation | Telkom University (sesuaikan) |
| Description | salin paragraf pertama README |
| Resource type | Software |
| License | MIT (kode). Tambahkan catatan CC BY 4.0 untuk CSV di deskripsi |
| Keywords | sponge hash, VSH, preimage, IoT, reproducibility |
| Related identifiers | nanti DOI artikel PeerJ, relasi `isSupplementTo` |
| Version | paper-v1 |
| Language | eng |

5. Publish. Zenodo memberi DOI `10.5281/zenodo.xxxxxxx`.
6. Salin DOI ke:
   - `CITATION.cff`
   - `README.md`
   - `DATA_AVAILABILITY.md`
   - naskah (kalimat Data Availability)
   - form online PeerJ

Versi konsep Zenodo memakai DOI versi (`.../zenodo.111`) dan DOI konsep (`.../zenodo.110`). Untuk paper, sitir **DOI versi `paper-v1`**.

---

## E. Jika GitHub-Zenodo tidak terhubung

Unggah zip manual:

1. Pack folder `spvsh320-repro` (tanpa `.git` juga boleh).
2. New upload di https://zenodo.org/uploads
3. Resource type Software
4. Isi metadata sama seperti tabel di atas
5. Publish → dapat DOI

Cara ini sah. Integrasi GitHub hanya memudahkan versi berikutnya.

---

## F. Form PeerJ

- **Data Availability**: tempel teks di `DATA_AVAILABILITY.md`.
- **Supplemental files**: jangan unggah stream 100 MB. Cukup DOI. File di situs PeerJ dibatasi 30 MB/file dan 50 MB total.
- Jangan unggah PDF tabel sebagai “raw data”. CSV di paket ini sudah machine-readable.

---

## G. Cek sebelum menekan Publish di Zenodo

- [ ] README menyebut tag `paper-v1`
- [ ] Tiga file `.cpp` dan sembilan file `run01` ada di zip
- [ ] Tidak ada password, kunci, atau path server pribadi
- [ ] Angka di `results/preimage_combined_summary.csv` sama dengan tabel naskah
- [ ] LICENSE ada
- [ ] Nama penulis sama dengan naskah
- [ ] Repo GitHub Public, atau staf jurnal bisa membuka tautan private

---

## H. Setelah DOI keluar, commit sekali lagi

```bash
# isi DOI di README.md, CITATION.cff, DATA_AVAILABILITY.md
git add README.md CITATION.cff DATA_AVAILABILITY.md
git commit -m "Record Zenodo DOI for paper-v1"
git push
```

Jangan mengganti tag `paper-v1` yang sudah diarsip. Jika DOI harus masuk ke snapshot yang sama, buat tag `paper-v1.1` dan versi Zenodo baru, lalu sitir yang baru di naskah revisi.
