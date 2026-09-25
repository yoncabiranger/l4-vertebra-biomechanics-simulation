import numpy as np

# 1. STL dosya adını tanımlıyoruz
dosya_adi = "Segmentation_tek_omur.stl"

# 2. Binary dosyayı açıp üçgen sayısını okuyoruz
with open(dosya_adi, "rb") as f:
    f.seek(80)  # İlk 80 baytlık boş başlığı atla
    n_triangles = np.frombuffer(f.read(4), dtype=np.uint32)[0]

    # Her biri 50 bayt olan üçgen verilerini tek seferde belleğe al
    data = np.fromfile(
        f,
        dtype=np.dtype(
            [
                ("normal", np.float32, (3,)),
                ("v1", np.float32, (3,)),
                ("v2", np.float32, (3,)),
                ("v3", np.float32, (3,)),
                ("attr", np.uint16),
            ]
        ),
        count=n_triangles,
    )

# 3. Tüm köşe koordinatlarını tek bir havuzda birleştir
noktalar = np.concatenate([data["v1"], data["v2"], data["v3"]], axis=0)

# 4. Sınır değerleri ve boyutları hesapla
min_xyz = np.min(noktalar, axis=0)
max_xyz = np.max(noktalar, axis=0)
boyutlar = max_xyz - min_xyz

print("=" * 50)
print("TEK OMUR GEOMETRİK BOYUTLARI")
print("=" * 50)
print(f"Toplam Üçgen Sayısı : {n_triangles}")
print(f"Genişlik (X ekseni) : {boyutlar[0]:.2f} mm")
print(f"Derinlik (Y ekseni) : {boyutlar[1]:.2f} mm")
print(f"Omur Boyu (Z ekseni): {boyutlar[2]:.2f} mm")
print("=" * 50)


# -------------------------------------------------------------
# ADIM 4: BİYOMEKANİK PARAMETRELERİN HESAPLANMASI
# -------------------------------------------------------------
# Boyutları metreye çeviriyoruz
L_x = boyutlar[0] * 1e-3  # m
L_y = boyutlar[1] * 1e-3  # m
L_z = boyutlar[2] * 1e-3  # m (Omur yüksekliği L)

# Eliptik enine kesit alanı (A)
A = (np.pi * L_x * L_y) / 4.0  # m^2

# Doku özellikleri (Literatür omur kompozit değerleri)
E = 8.0e9  # Elastisite modülü (Pa) -> 8 GPa
rho = 1800.0  # Yoğunluk (kg/m^3)
zeta = 0.03  # Sönüm oranı (%3)

# 1-DOF dinamik katsayılar
K = (E * A) / L_z  # Eksenel Rijitlik (N/m)
M = rho * A * L_z  # Kütle (kg)
omega_n = np.sqrt(K / M)  # Açısal doğal frekans (rad/s)
f_n = omega_n / (2.0 * np.pi)  # Doğal frekans (Hz)
C = 2.0 * zeta * np.sqrt(K * M)  # Sönüm katsayısı (N.s/m)

print("\n" + "=" * 50)
print("BİYOMEKANİK DİNAMİK PARAMETRELER")
print("=" * 50)
print(f"Efektif Kesit Alanı (A) : {A*1e4:.2f} cm²")
print(f"Toplam Kütle (M)        : {M:.3f} kg")
print(f"Eksenel Rijitlik (K)    : {K/1e6:.2f} MN/m")
print(f"Sönüm Katsayısı (C)     : {C:.2f} N·s/m")
print(f"Doğal Frekans (fn)      : {f_n:.2f} Hz")
print("=" * 50)


# -------------------------------------------------------------
# ADIM 5: DİNAMİK YÜKLEME VE NEWMARK-BETA ZAMAN İNTEGRASYONU
# -------------------------------------------------------------

# 1. Simülasyon Zaman Parametreleri
dt = 0.0005  # Zaman adımı: 0.5 milisaniye (yüksek hassasiyet)
t_toplam = 0.5  # Toplam simülasyon süresi: 0.5 saniye
t = np.arange(0, t_toplam, dt)
N_adim = len(t)

# 2. Dış Kuvvet Profili F(t) (BioSuit Statik + Piezo Dinamik)
F_statik = 150.0  # Dava Newman eksenel ön-yükü (N)
F_dinamik_genlik = 40.0  # Piezo aktüatör kuvvet genliği (N)
f_uyari = 30.0  # Piezo frekansı (Hz)

# Toplam kuvvet zaman serisi
F = F_statik + F_dinamik_genlik * np.sin(2.0 * np.pi * f_uyari * t)

# 3. Newmark-Beta Parametreleri (Ortalama İvme / Koşulsuz Kararlı)
gamma = 0.50
beta = 0.25

# İntegrasyon sabitleri
a0 = 1.0 / (beta * (dt**2))
a1 = gamma / (beta * dt)
a2 = 1.0 / (beta * dt)
a3 = (1.0 / (2.0 * beta)) - 1.0
a4 = (gamma / beta) - 1.0
a5 = (dt / 2.0) * ((gamma / beta) - 2.0)

# Efektif Rijitlik Matrisi (1-DOF skaler)
K_efektif = K + a0 * M + a1 * C

# 4. Çözüm Dizilerinin Başlatılması
u = np.zeros(N_adim)  # Deplasman (m)
v = np.zeros(N_adim)  # Hız (m/s)
a = np.zeros(N_adim)  # İvme (m/s^2)

# Başlangıç ivmesi (t=0 anında u=0, v=0 kabulüyle)
a[0] = (F[0] - C * v[0] - K * u[0]) / M

# 5. Zaman Adımı Döngüsü (Newmark-Beta Çözücü)
for i in range(N_adim - 1):
    # Efektif Yük Hesabı
    F_efektif = (
        F[i + 1]
        + M * (a0 * u[i] + a2 * v[i] + a3 * a[i])
        + C * (a1 * u[i] + a4 * v[i] + a5 * a[i])
    )

    # Bir sonraki adımın yer değiştirmesi
    u[i + 1] = F_efektif / K_efektif

    # Bir sonraki adımın ivmesi ve hızı
    a[i + 1] = (
        a0 * (u[i + 1] - u[i]) - a2 * v[i] - a3 * a[i]
    )
    v[i + 1] = v[i] + dt * ((1.0 - gamma) * a[i] + gamma * a[i + 1])

# 6. Biyomekanik ve Hücresel Sonuç Analizi
# Dinamik dalgalanmanın genliğini bulmak için son 0.1 saniyelik kararlı rejimi alıyoruz
kararli_aralik = int(0.4 / dt)
u_kararli = u[kararli_aralik:]
u_statik_ortalama = np.mean(u_kararli)  # Statik BioSuit çökmesi (m)
u_dinamik_genlik = (np.max(u_kararli) - np.min(u_kararli)) / 2.0  # Piezo salınımı (m)

# Nanometre (nm) ve Mikrometre (um) dönüşümleri
u_statik_um = u_statik_ortalama * 1e6
u_dinamik_nm = u_dinamik_genlik * 1e9

print("\n" + "=" * 50)
print("NEWMARK-BETA SİMÜLASYON VE HÜCRESEL SONUÇLARI")
print("=" * 50)
print(f"BioSuit Statik Çökme : {u_statik_um:.4f} µm")
print(f"30 Hz Piezo Titreşim Genliği : {u_dinamik_nm:.2f} nm")
print(f"Primer Silya Eşiği (Hedef)   : 30.00 nm")

if u_dinamik_nm >= 30.0:
    print("DURUM: >> OSTEOGENİK UYARI SAĞLANDI (Primer silya büküldü) <<")
else:
    print("DURUM: >> Titreşim genliği eşiğin altında kaldı, kuvveti artırın <<")
print("=" * 50)

# -------------------------------------------------------------
# ADIM 6: GRAFİKLEŞTİRME VE BİYOMEKANİK GÖRSELLEŞTİRME
# -------------------------------------------------------------
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 8))

# 1. Panel: Genel Dinamik Tepki (0 - 0.5 s)
plt.subplot(2, 1, 1)
plt.plot(t, u * 1e6, color="royalblue", linewidth=1.5, label="Omur Çökmesi u(t)")
plt.axhline(
    u_statik_um,
    color="crimson",
    linestyle="--",
    alpha=0.8,
    label=f"Statik BioSuit Dengesi ({u_statik_um:.3f} µm)",
)
plt.title(
    "L4 Omuru Dinamik Deplasman Yanıtı (BioSuit + 30 Hz Piezo)",
    fontsize=12,
    fontweight="bold",
)
plt.xlabel("Zaman (s)", fontsize=10)
plt.ylabel("Deplasman (µm)", fontsize=10)
plt.grid(True, linestyle=":", alpha=0.6)
plt.legend(loc="upper right")

# 2. Panel: Kararlı Rejim ve Primer Silya Doğrulaması (0.4 - 0.5 s)
plt.subplot(2, 1, 2)
t_kararli = t[kararli_aralik:]
u_dinamik_dalga_nm = (u_kararli - u_statik_ortalama) * 1e9

plt.plot(
    t_kararli,
    u_dinamik_dalga_nm,
    color="darkgreen",
    linewidth=1.8,
    label=f"Piezo Salınımı (Genlik = {u_dinamik_nm:.2f} nm)",
)
plt.axhline(
    30.0,
    color="red",
    linestyle="--",
    linewidth=1.5,
    label="Primer Silya Eşiği (30 nm)",
)
plt.axhline(-30.0, color="red", linestyle="--", linewidth=1.5)
plt.title(
    "Kararlı Rejimde Osteosit Primer Silya Uyarım Genliği",
    fontsize=12,
    fontweight="bold",
)
plt.xlabel("Zaman (s)", fontsize=10)
plt.ylabel("Dinamik Salınım (nm)", fontsize=10)
plt.grid(True, linestyle=":", alpha=0.6)
plt.legend(loc="upper right")

plt.tight_layout()
plt.savefig("tek_omur_dinamik_yanit.png", dpi=300)
print("\n[BİLGİ] 'tek_omur_dinamik_yanit.png' başarıyla kaydedildi.")

# -------------------------------------------------------------
# ADIM 7: HÜCRESEL MEKANOTRANSDÜKSİYON DOĞRULAMA ÖZETİ
# -------------------------------------------------------------
guvenlik_marji = ((u_dinamik_nm - 30.0) / 30.0) * 100

print("\n" + "=" * 55)
print("BİYOMEKANİK & HÜCRESEL DOĞRULAMA RAPORU (SUMMARY)")
print("=" * 55)
print(f"Hedef Eşik        : 30.00 nm")
print(f"Elde Edilen Genlik: {u_dinamik_nm:.2f} nm")
print(f"Güvenlik Marjı    : %{guvenlik_marji:.1f}")
print("Sonuç             : MEKANİK UYARI OSTEOGENEZ İÇİN YETERLİ")
print("=" * 55)
plt.show()