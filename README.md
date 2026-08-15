# registry-distress-monitor

Izlog zatvorenog sustava koji vremenski nadzire hrvatske javne registre i trazi
rani znak da subjekt ili nekretnina ulazi u nevolju. Objavljen je jedan modul u
cijelosti: konformna kalibracija cjenovnog raspona.

Pisano kao pitanja i odgovori, jer su to pitanja koja se zaista postavljaju.

---

**Sto tocno nadzire?**

Blokade racuna, stecajne objave, javne drazbe i promjene u sudskom registru,
povezane po identifikatoru pravne osobe i katastarskoj cestici. Uzi izdanak,
cjevovod prikupljanja s ponovljivim mjerenjem modela, objavljen je zasebno kao
[croatian-registry-pipeline](https://github.com/hm53-byte/croatian-registry-pipeline).

**Koliko je toga stvarno napisano?**

Mjereno 15. 8. 2026. na Windowsu 11, Python 3.13.1: 2089 testova prolazi, 2
preskocena, nijedan ne pada, u 193 sekunde. Kod ima 41 370 redaka u 176 modula,
testovi 35 550 redaka u 158 datoteka. Migracija sheme je 62.

Odnos koda i testova je oko 1 prema 0,86. Sustav donosi tvrdnje o tudjem
poslovanju na temelju javnih zapisa, pa je cijena tihe pogreske veca od cijene
testa.

**Mogu li to provjeriti?**

Ne za brojke iznad; one su iz zatvorenog repozitorija. Mozes provjeriti
objavljeni modul:

```bash
git clone https://github.com/hm53-byte/registry-distress-monitor
cd registry-distress-monitor
pip install pytest && python -m pytest testovi -q      # 21 test
```

**Je li nesto od zamisljenog propalo?**

Tri stvari, i sve tri stoje zapisane.

Prva. Sredisnja zamisao rane inacice bila je da presusenje aktivnosti prethodi
nevolji. Mjereno na povijesnim podacima, prognosticka vrijednost je bila lift
0,12 puta uz p = 1,0 na horizontu od 30 do 45 dana. To nije slab rezultat nego
nikakav: signal je losiji od slucajnog. Uzrok je strukturan. Praznina u zapisu
mnogo cesce znaci da izvor nije objavio nego da se nista nije dogodilo, pa se
mjeri pokrivenost izvora, a ne stanje subjekta. Klasa signala utemeljenih na
odsutnosti je napustena.

Druga. Odziv je dugo bio nizak i pretpostavka je bila da model nije dovoljno
dobar. Mjereno odvojeno, najveci gubitak nastaje **prije** modela: dio objava
uopce ne udje u bazu. Model koji radi nad polovicom dogadjaja ne popravlja se
boljim modelom.

Treca. Predlozen je Bloomov filtar ispred provjere postojanja zapisa, izmjereno
ubrzanje 2,8 puta. Provjerom se pokazalo da vrijedi samo za promasaje, a da nad
stupcem s jedinstvenim indeksom SQLite vec radi dovoljno brzo. Filtar je ostao,
ali iskljucen po zadanom, s napisanim uvjetom pod kojim se isplati.

**Sto radi objavljeni modul?**

Racuna raspon cijene nad malim uzorkom, gdje obicna procjena laze u oba smjera.
Tri primitiva pokrivaju tri rezima velicine uzorka.

Za n >= 30 ide Mondrian konformna predikcija, bez pretpostavke o distribuciji,
uz pokrivenost barem 1 - alpha:

```python
scores = sorted(abs(x - med) for x in amounts)
rank = math.ceil((n + 1) * (1 - alpha))   # Vovkova korekcija
rank = min(rank, n)                        # bez ovoga indeks izadje iz polja
q = scores[rank - 1]
```

Za 3 <= n < 30 dolazi James-Stein sazimanje, ali samo ako KS test propusti. Ako
se lokalna distribucija razlikuje od globalne, globalni prosjek nije
reprezentativan i sazimanje bi uvelo pristranost umjesto da je smanji, pa se
pada na Wilsonov interval. Razlog pada se zapisuje, da se poslije moze
razlikovati "nije bilo globalnog uzorka" od "globalni uzorak nije bio
usporediv".

Ispod n = 3 nema pojasa. Vraca se `no_data`, ne uzak pojas s laznom sigurnoscu.

**Ima li modul poznatih kvarova?**

Ima jedan i objavljen je. Testovi pisani za ovaj izlog pokazali su da put
oznacen kao `james_stein` zapravo **ne sazima**: ulazna tocka predaje formuli
jednu grupu, a formula trazi barem tri, pa vraca ulaz nepromijenjen i faktor
1,0. Vraceni medijan jednak je lokalnom medijanu.

Modul je objavljen kakav jest, uz test `test_sazimanje_u_dual_path_ne_saziva`
koji to pribija, da promjena bude vidljiva a ne tiha. Mondrianov put, koji je
glavni, tim nalazom nije dotaknut. Popravak je poznat i nije ugradjen:
sazimanje ima smisla tek kad se proslijede sve katastarske opcine odjednom kao
grupe.

**Kako se rjesava pravni okvir?**

Kodom, ne dogovorom. Ruta pretrage po imenu osobe odbija se prije nego zahtjev
nastane, i provjera se izvodi dvaput, nad putanjom i nad sastavljenim naslovom.
Ritam prema posluzitelju ima dvije neovisne brane i brojac zapreka koji prezivi
gasenje procesa. Izvor ciji uvjeti nisu procitani vodi se kao neprovjeren, ne
kao dopusten. Osobni podaci se rastavljaju u memoriji, ali za njih u shemi nema
stupca, pa na disk ne izlaze.

**Sto nije objavljeno?**

Kod za prikupljanje, shema baze, nijedan prikupljeni zapis i nijedan naziv
subjekta.

---

Apache-2.0, [LICENSE](LICENSE).
