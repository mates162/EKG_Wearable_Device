# ESP32 12‑svodové EKG — Real-Time WiFi Plotter (v11)

Desktopová aplikace pro živé zobrazení a záznam 12‑svodového EKG z ESP32 přes WiFi. Přijímá binární data z ESP32 (ADS1298, 8 kanálů @ 500 SPS), dopočítává odvozené svody (III, aVR, aVL, aVF) a zobrazuje plné 12‑svodové EKG s filtry, záznamem do CSV a odhadem tepové frekvence.

---

## Požadavky a spuštění

- **Python:** 3.x  
- **Závislosti:** `numpy`, `pyqtgraph`, `PyQt5`, `scipy`  
- **Instalace:**  
  `pip install -r requirements.txt`  
  nebo:  
  `pip install pyqtgraph PyQt5 numpy scipy`

**Spuštění:**

1. Připojte PC k WiFi přístupovému bodu ESP32 (SSID např. `ESP32_ECG`, heslo dle firmware).
2. Spusťte:  
   `python esp_save_plot_V11.py`
3. Aplikace se sama připojí na `192.168.4.1`, port **3333** (data), **3334** (příkazy) a spustí stream.

---

## Přehled funkcí

### 1. Živý příjem a zobrazení EKG

- **Připojení:** TCP klient na ESP32 (výchozí `192.168.4.1:3333`).
- **Protokol:** Binární rámce `[0xA5][0x5A][počet][vzorky…]`, každý vzorek = 4 B časová známka (µs, LE) + 27 B surová data ADS1298 (3 B status + 8×3 B kanály, 24‑bit signed big‑endian).
- **Validace:** Pouze vzorky s platným statusem ADS1298 (horní nibble `0xC`) se zpracují; neplatné se přeskočí.
- **12 svodů:** Z 8 surových kanálů se dopočítají III, aVR, aVL, aVF; zobrazení v pořadí **I, II, III, aVR, aVL, aVF, V1–V6**.
- **Ring buffer:** Posledních **10 s** dat (5000 vzorků při 500 Hz) v paměti pro plynulé vykreslování.
- **Plynulý scroll:** Osa X plynule sleduje „nejnovější“ čas (nastavitelné alpha).
- **OpenGL:** Vykreslování s OpenGL pro plynulý obraz (~60 FPS).

---

### 2. Klinické filtry EKG

- **Pásmový filtr (AAMI EC11):** 0,5–40 Hz, Butterworth 4. řádu (baseline wander + vysokofrekvenční šum).
- **Notch filtr:** 50 Hz nebo 60 Hz (síťové rušení), volba v comboboxu.
- **Okraje:** Delší padding na začátku/konci signálu (`FILTER_EDGE_PAD_SAMPLES`) pro menší deformaci na kraji dat.
- **Režim Y při zapnutých filtrech:**  
  - Výchozí fixní rozsah **−1 až +1 mV**.  
  - Zoom kolečkem jen **symetricky kolem Y = 0** (autoscale neodjíždí od nuly).  
  - Tlačítko **„Autoscale Y“** vrátí rozsah −1..1 mV (pouze při zapnutých filtrech).
- **Režim Y při vypnutých filtrech:** Plný autoscale osy Y, zoom/pan bez omezení.

Filtry vyžadují **scipy**; bez něj se data zobrazí nefiltrovaná.

---

### 3. Osa Y a zoom

- **Při zapnutých filtrech:**  
  - Fixní rozsah −1..1 mV (obnovuje se periodicky, pokud uživatel nezoomoval).  
  - Při zoomu kolečkem nebo tažením se vždy vynutí symetrický rozsah kolem 0 (`[-span, +span]`), min. span 0,02 mV, max. 50 mV.
- **Při vypnutých filtrech:** Automatické škálování Y (autoscale), volný zoom a pan.
- **Synchronizace X:** Všechny kanály sdílí stejný rozsah osy X (zoom/pan X společný).

---

### 4. Záznam do CSV (REC)

- **Tlačítko REC / STOP:** Spustí/zastaví záznam.
- **Obsah souboru:** Jedna řádka na vzorek: `timestamp_us`, 12 sloupců v mV (I, II, III, aVR, aVL, aVF, V1–V6), `gain`.
- **Kde se ukládá:** Do složky `logs` vedle skriptu (nebo zvolená cesta); výchozí název `ecg_data_YYYYMMDD_HHMMSS.csv`.
- **Vzorkovací frekvence:** Plných **500 SPS** (žádné podvzorkování); záznam běží v příjmovém vlákně nezávisle na vykreslování.

---

### 5. Načtení a prohlížení CSV

- **Načíst CSV:** Tlačítko **„Načíst CSV“** otevře dialog; očekávaný formát: sloupec `timestamp_us` a 12 kanálů v mV (stejné pořadí jako při záznamu).
- **Režim zobrazení:** Combobox **„Živý signál“** / **„Načtený CSV“**. V režimu CSV se zobrazuje načtený soubor.
- **Posun v čase:** Slider posouvá **10s okno** v rámci načteného záznamu; pod sliderem je text s aktuálním časem a počtem vzorků.
- **Filtry a Y:** Na načtená data se stejně aplikují klinické filtry (zapnuto/vypnuto) a chování osy Y (fixní −1..1 mV při zapnutých filtrech, symetrický zoom kolem 0).

---

### 6. Odhad tepové frekvence (HR)

- **Metoda:** Detekce R‑peaků z energie derivace signálu (prahování), odhad BPM z intervalu RR.
- **Kanály:** Automatický výběr nejlepšího kanálu z končetinových svodů (I, II) pro stabilitu odhadu (kritérium kvality z variability RR).
- **Výstup:** V status řádku např. `HR: 72 BPM (II)`; při CSV režimu stejný odhad pro aktuální výřez.
- **Aktualizace:** Přepočet cca 1× za sekundu (throttling).

---

### 7. Dvojklik – maximalizace kanálu

- **Dvojklik** na jeden z 12 subplotů: daný kanál se **maximalizuje** (zbytek skryt).
- **Dvojklik znovu** na maximalizovaný kanál: návrat na **všechny kanály** (normální výška, synchronizovaná X).

---

### 8. Pozastavení vykreslování

- Tlačítko **„Pozastavit vykreslování“** / **„Spustit vykreslování“**: zastaví/spustí timer, který aktualizuje křivky a osy. Data dál přicházejí do ring bufferu; po spuštění se zobrazí aktuální stav.

---

### 9. Dálkové ovládání zesílení (GAIN) na ESP32

- **Combobox GAIN:** Volby 1, 2, 3, 4, 6, 8, 12 (platné hodnoty PGA ADS1298).
- **Příkaz:** Po změně se na port **3334** odešle řádek `GAIN:<hodnota>\n`; očekává se odpověď typu `GAIN_OK:<potvrzená_hodnota>`.
- **Důsledky:**  
  - Aktualizuje se převod LSB → mV pro nové vzorky.  
  - Ring buffer se po potvrzení GAIN vyprázdní (staré vzorky byly v jiném zesílení).  
  - V status řádku se zobrazí např. „GAIN x6 potvrzeno ✓“ nebo chybová hláška.

---

### 10. Status řádek

- **Živý režim:** Stav připojení (🟢/🔴), počet přijatých vzorků, počet rámců, vzorky/s, zaplnění bufferu, latence (ms), jitter časových známek (µs), aktuální GAIN, HR (BPM + kanál).
- **CSV režim:** Název souboru, časové okno, počet vzorků v okně, HR pro výřez.

---

### 11. Uzavření aplikace

- Při zavření okna se zastaví příjmové vlákno a při běžícím záznamu se CSV soubor korektně uzavře (flush + close).

---

## Konfigurace (konstanty v kódu)

| Konstanta | Význam | Výchozí |
|-----------|--------|---------|
| `ESP_HOST` | IP ESP32 (AP) | `192.168.4.1` |
| `TCP_PORT` | Port datového streamu | `3333` |
| `TCP_PORT_CMD` | Port příkazů (GAIN) | `3334` |
| `SAMPLE_RATE` | Vzorkovací frekvence [Hz] | `500.0` |
| `WINDOW_SEC` | Délka zobrazeného okna [s] | `10.0` |
| `PLOT_INTERVAL_MS` | Interval obnovy grafu [ms] | `12` |
| `Y_VIEW_MIN_MV`, `Y_VIEW_MAX_MV` | Fixní rozsah Y při filtrech [mV] | `-1.0`, `1.0` |
| `FILTER_BANDPASS_LOW/HIGH` | Pásmo filtru [Hz] | `0.5`, `40.0` |
| `FILTER_EDGE_PAD_SAMPLES` | Délka paddingu na okrajích [vzorky] | `1.2 * SAMPLE_RATE` |

---

## Struktura projektu (relevantní soubory)

- **`esp_save_plot_V11.py`** – hlavní skript plotru (příjem, filtry, GUI, záznam, načtení CSV, HR).
- **`requirements.txt`** – závislosti Pythonu.
- **`src/main.cpp`** – firmware ESP32 (protokol a rozhraní musí odpovídat v11).
- **`logs/`** – výchozí adresář pro ukládání CSV záznamů.

---

## Shrnutí funkcí v bodech

1. Živý příjem 12‑svodového EKG přes WiFi (TCP, binární protokol v11).  
2. Validace statusu ADS1298, výpočet odvozených svodů a převod na mV.  
3. Klinické filtry: pásmový 0,5–40 Hz + notch 50/60 Hz s menší deformací na okrajích.  
4. Osa Y: při filtrech fixní −1..1 mV a zoom jen symetricky kolem Y=0; bez filtrů plný autoscale.  
5. Záznam do CSV v plném 500 SPS.  
6. Načtení CSV a prohlížení s posuvníkem a stejnými filtry.  
7. Odhad tepové frekvence s auto‑výběrem kanálu.  
8. Dvojklik pro maximalizaci jednoho kanálu.  
9. Pozastavení/spuštění vykreslování.  
10. Dálková změna GAIN na ESP32 přes příkazový port.  
11. Status řádek s připojením, statistikami a HR.
