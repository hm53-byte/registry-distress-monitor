# registry-distress-monitor

Izlog jednog privatnog sustava: vremenski nadzor hrvatskih javnih registara
koji trazi rani znak da poslovni subjekt ili nekretnina ulazi u nevolju.
Blokade racuna, stecajne objave, javne drazbe, promjene u sudskom registru,
sve povezano po identifikatoru pravne osobe i katastarskoj cestici.

Ovdje nije cijeli sustav. Ovdje su brojke, tri nalaza koja su promijenila
smjer projekta, i **jedan cijeli modul s testovima**: konformna kalibracija
cjenovnog raspona.

Uzi izdanac ovog rada, cjevovod prikupljanja s ponovljivim mjerenjem modela,
objavljen je zasebno:
[croatian-registry-pipeline](https://github.com/hm53-byte/croatian-registry-pipeline).

---

## Izmjereno, ne prepisano

Pokrenuto nad izvornim repozitorijem 15. 8. 2026., Windows 11, Python 3.13.1:

| Mjera | Vrijednost |
|---|---|
| Testova | **2089 prolazi, 2 preskocena, 0 padova** |
| Trajanje skupa | 193 s |
| Redaka Pythona u `src/` | 41 370 u 176 modula |
| Redaka testova | 35 550 u 158 datoteka |
| Migracija sheme | 62 |

Odnos koda i testova je otprilike 1 prema 0,86. To nije slucajno: sustav
donosi tvrdnje o tudjem poslovanju na temelju javnih zapisa, pa je cijena
tihe pogreske veca od cijene testa.

---

## Tri nalaza koja su promijenila smjer

Ovo je dio koji izlog obicno preskoci. Sva tri su negativna i sva tri su
zadrzana.

### 1. Odsutnost signala nije signal

Sredisnja zamisao rane inacice bila je da presusenje aktivnosti, dakle
razdoblje bez ijedne objave, prethodi nevolji. Mjereno na povijesnim
podacima, prognosticka vrijednost je bila **lift 0,12 puta uz p = 1,0** na
horizontu od trideset do cetrdeset pet dana.

To nije slab rezultat nego nikakav: signal je losiji od slucajnog. Uzrok je
strukturan, ne statisticki. Praznina u zapisu mnogo cesce znaci da izvor nije
objavio nego da se nista nije dogodilo, pa se mjeri pokrivenost izvora, a ne
stanje subjekta.

Posljedica: klasa signala utemeljenih na odsutnosti je napustena, a mjera
pokrivenosti izvora uvedena kao zaseban pokazatelj.

### 2. Usko grlo je bio unos, ne model

Odziv sustava dugo je bio nizak i pretpostavka je bila da model nije dovoljno
dobar. Kad je mjereno odvojeno, pokazalo se da najveci gubitak nastaje prije
modela: dio objava uopce ne udje u bazu, zbog oblika dokumenta ili izostanka
poveznice s pravnom osobom.

Model koji radi nad polovicom dogadjaja ne moze se popraviti boljim modelom.
Trud se preselio na unos.

### 3. Postojeci indeks vec je bio brz

Predlozen je Bloomov filtar ispred provjere postojanja zapisa i izmjereno je
ubrzanje od 2,8 puta. Provjerom se pokazalo da ubrzanje vrijedi samo
za promasaje, a da nad stupcem s jedinstvenim indeksom SQLite vec radi
dovoljno brzo za ovu kolicinu podataka.

Filtar je ostao u kodu, ali iskljucen po zadanom, s napisanim uvjetom pod
kojim se isplati (mrezni poziv iza provjere, ne lokalni upit).

---

## Objavljeni modul: konformna kalibracija cjenovnog raspona

[`primjer/conformal_calibration.py`](primjer/conformal_calibration.py), 254
retka, bez ijedne obavezne ovisnosti. `scipy` se koristi ako postoji, inace
radi vlastita izvedba KS testa.

Problem: za raspon cijene nad malim uzorkom obicna procjena laze u oba
smjera. Tri primitiva pokrivaju tri rezima velicine uzorka.

**Za n ≥ 30, Mondrian konformna predikcija.** Bez pretpostavke o distribuciji,
uz pokrivenost barem 1 − alpha:

```python
scores = sorted(abs(x - med) for x in amounts)
rank = math.ceil((n + 1) * (1 - alpha))   # Vovkova korekcija
rank = min(rank, n)                        # bez ovoga indeks izadje iz polja
q = scores[rank - 1]
```

Ogranicenje ranga na `n` nije kozmeticko: za malen alpha izraz
`(n+1)(1-alpha)` prelazi `n`. Test to pribija.

**Za 3 ≤ n < 30, James-Stein sazimanje, ali samo ako brana propusti.** Prije
sazimanja se KS testom provjerava razlikuje li se lokalna distribucija od
globalne. Ako se razlikuje, globalni prosjek nije reprezentativan i sazimanje
bi uvelo pristranost umjesto da je smanji, pa se pada na Wilsonov interval.
Razlog pada zapisuje se u rezultat (`ks_gated`), da se poslije moze
razlikovati "nije bilo globalnog uzorka" od "globalni uzorak nije bio
usporediv".

**Ispod n = 3 nema pojasa.** Vraca se `no_data`, ne uzak pojas s laznom
sigurnoscu.

### Sto je mjerenje pokazalo o samom modulu

Testovi pisani za ovaj izlog otkrili su da **put oznacen kao `james_stein`
zapravo ne sazima**. Ulazna tocka predaje formuli jednu jedinu grupu, a
formula trazi barem tri, pa vraca ulaz nepromijenjen i faktor 1,0. Vraceni
medijan jednak je lokalnom medijanu; globalni prosjek ne utjece ni na sto.

Modul je objavljen kakav jest, uz test koji to ponasanje pribija
(`test_sazimanje_u_dual_path_ne_saziva`), da promjena bude vidljiva, a ne
tiha. Mondrianov put, koji je glavni, tim nalazom nije dotaknut.

Popravak je poznat i nije ugradjen ovdje: sazimanje ima smisla tek kad se
proslijede sve katastarske opcine odjednom kao grupe, a ne jedna po jedna.

---

## Pravni okvir kao dio arhitekture

Sustav dohvaca javne registre, pa granice nisu stvar dogovora nego koda:

- Zabrana pojedine rute na razini putanje, ne samo domacina. Ruta pretrage po
  imenu osobe odbija se prije nego zahtjev nastane, i provjera se izvodi
  dvaput, nad putanjom i nad sastavljenim naslovom.
- Ritam prema posluzitelju s dvije neovisne brane i brojacem zapreka koji
  prezivi gasenje procesa.
- Izvor ciji uvjeti nisu procitani vodi se kao neprovjeren, a ne kao dopusten.
  Sutnja nije dozvola.
- Osobni podaci se rastavljaju u memoriji, ali za njih u shemi nema stupca,
  pa na disk ne izlaze.

---

## Sto ovdje nije

Nije objavljen kod za prikupljanje, shema baze, nijedan prikupljeni zapis ni
ijedan naziv subjekta. Sve brojke izmjerene su nad zatvorenim repozitorijem.

---

## Pokretanje ovog repozitorija

```
python -m pytest testovi -q
```

21 test, samo standardna biblioteka i `pytest`.

---

## Licenca

Apache-2.0, vidi [LICENSE](LICENSE).
