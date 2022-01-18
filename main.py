import numpy as np

import solsystem as rita


def getacc(m, p, G, N):
    a = np.zeros((N, 3), dtype='g')
    for i in range(N):
        for j in range(N):
            dx = np.longfloat(p[j, 0] - p[i, 0])
            dy = np.longfloat(p[j, 1] - p[i, 1])
            dz = np.longfloat(p[j, 2] - p[i, 2])
            inv_r3 = np.longfloat((dx ** 2) + (dy ** 2) + (dz ** 2) +
                                  (0.01 ** 2)) ** -1.5
            a[i, 0] += G * (dx * inv_r3) * m[j]
            a[i, 1] += G * (dy * inv_r3) * m[j]
            a[i, 2] += G * (dz * inv_r3) * m[j]
    return a


def getKE(m, v, N):
    vt = np.zeros((N, 1), dtype='g')
    KE = np.zeros((N, 1), dtype='g')
    KEt = np.longfloat(0.0)
    for i in range(N):
        vt[i] = (v[i, 0] ** 2 + v[i, 1] ** 2 + v[i, 2] ** 2) ** 0.5
        KE[i] = 0.5 * (m[i] * vt[i] ** 2)
        KEt += KE[i]
    return KEt[0]


def getPE(m, p, G, N):
    PE = np.longfloat(0.0)
    inv_r = np.zeros((N, 1), dtype='g')
    for i in range(N):
        for j in range(N):
            if i != j:
                dx = np.longfloat(p[j, 0] - p[i, 0])
                dy = np.longfloat(p[j, 1] - p[i, 1])
                dz = np.longfloat(p[j, 2] - p[i, 2])
                inv_r[i] = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                inv_r[i] = 1.0 / inv_r[i]
                PE += G * (-(m[i] * m[j]) * inv_r[i])
    return 0.5 * PE[0]


print("Detta är en tyngdkrafts simulation skapad av Arvid och Isak")
N = 6
# antalet kroppar i simulationen
m = np.zeros((N, 1), dtype='g')
# 1-d array som innehåller massan av kropparna efter index N i enheten kilogram
m[0] = 1.9885 * 10 ** 30
# massan av solen 1,988,500*10**24kg 1.9885*10**30
m[1] = 3.3011 * 10 ** 23
# merkurius
m[2] = 4.8675 * 10 ** 24
# venus
m[3] = 5.9724 * 10 ** 24
# massan av jorden 5.9724 * 10**24kg
m[4] = 7.346 * 10 ** 22
# massan av månen 0.07346*10**24kg 7.346*10**22kg
m[5] = 6.4171 * 10 ** 23
# mars
p = np.zeros((N, 3), dtype='g')
# En 2-d array med formen N*[0,0,0] och enheten meter,
# p[0] är positionen i x-led, p[1] är positionen i y-led och p[2] är positionen i z-led.
p[0] = [0, 0, 0]
# pos. solen
p[1] = [4.6002 * 10 ** 10, 0, 0]
p[2] = [1.07476 * 10 ** 11, 0, 0]
p[3] = [1.47092 * 10 ** 11, 0, 0]
# Om jorden börjar vid avståndet [perihelion=x, 0, 0],
# så måste dess hastighet vara riktad mot y med vinkeln i mot z på x-y planet.
# pos. jorden 147.092*10**6km 1.47092*10**11m
p[4] = [1.47092 * 10 ** 11, 3.633 * 10 ** 8, 0]  # pos. månen y=0.3633*10**6km från jorden 3.633*10**8m
p[5] = [2.06617 * 10 ** 11, 0, 0]
v = np.zeros((N, 3), dtype='g')
# En 2-d array med formen N*[0,0,0] och enheten meter per sekund,
# v[0] är hastigheten i x-led, v[1] är positionen i y-led och v[2] är positionen i z-led.
v[0] = [0, 0, 0]
v[1] = [0, 5.85397445423 * 10 ** 4, 7.19296245783 * 10 ** 3]  # i=7.005 v=5.898*10^4 m/s
v[2] = [0, 3.51981186852 * 10 ** 4, 2.08807112494 * 10 ** 3]  # i=3.395 v=3.526*10^4 m/s
v[3] = [0, 3.029 * 10 ** 4, 0]
# hastighet y jorden 30.29km/s 3.029*10^4 m/s i = 0 (deg).
v[4] = [-1.07764055794 * 10 ** 3, 3.029 * 10 ** 4, 9.7030035964 * 10 ** 1]
# hastighet v-xy månen -1.082km/s -1.082*10^3 m/s i = 5.145 (deg). x*x + z*z = v*v z = v * sin(i)
v[5] = [0, 2.6486172436 * 10 ** 4, 8.55961268228 * 10 ** 2]  # i=1.851 v=2.65*10^4 m/s
G = np.longfloat(6.6743 * 10 ** -11)

t = 0.0
tslut = 31536000
# tslut = 63072000
dt = 86400
# dt = 3600
Nt = int(
    np.ceil(tslut / dt)
)
a = getacc(m, p, G, N)
KE = getKE(m, v, N)
PE = getPE(m, p, G, N)

p_save = np.zeros((N, 3, Nt + 1), dtype=[('pos', 'g'), ('tid', 'l')])
p_save[:, :, 0]['pos'] = p
KE_save = np.zeros(Nt + 1, dtype='g')
KE_save[0] = KE
PE_save = np.zeros(Nt + 1, dtype='g')
PE_save[0] = PE

for it in range(Nt):
    for i in range(N):
        v[i, 0] += a[i, 0] * (dt / 2.0)
        v[i, 1] += a[i, 1] * (dt / 2.0)
        v[i, 2] += a[i, 2] * (dt / 2.0)
        p[i, 0] += v[i, 0] * dt
        p[i, 1] += v[i, 1] * dt
        p[i, 2] += v[i, 2] * dt
    a = getacc(m, p, G, N)
    for i in range(N):
        v[i, 0] += a[i, 0] * (dt / 2.0)
        v[i, 1] += a[i, 1] * (dt / 2.0)
        v[i, 2] += a[i, 2] * (dt / 2.0)
    t += dt
    KE = getKE(m, v, N)
    PE = getPE(m, p, G, N)
    p_save[:, :, it + 1]['pos'] = p
    p_save[:, :, it + 1]['tid'] = dt * (it + 1)
    KE_save[it + 1] = KE
    PE_save[it + 1] = PE
# TE = KE_save + PE_save  #Räknar ut total energin och tar fram medelvärdet
# print(TE[0])
# print(np.mean(TE, 0))
# print(TE)
# px = np.zeros((N,Nt+1))
# py = np.zeros((N,Nt+1))
# pz = np.zeros((N,Nt+1))
# print(np.shape(px))
# print(np.shape(p_save))
# px = p_save[:,0,:]
# print(np.shape(px))
# py = p_save[:,1,:]
# pz = p_save[:,2,:]
dt_maxmin = np.dtype([
    ('pos', 'g'), ('tid', 'l')
])  # tid blir i enheten dt, då den bara sparar informationen varje tidssteg.
p_max = np.zeros((N, 3), dtype=dt_maxmin)
# print(np.shape(px_max))
# print(px_max)
p_min = np.zeros((N, 3), dtype=dt_maxmin)
# print(px)
# print(py)
# print(pz)

# print((np.amax(p_save[2,2,:]['pos'],0)))
# p_max[2,2]['pos']=(np.amax(p_save[2,2,:]['pos'],0))
# print(p_max[2,2]['pos'])
# print(p_save[2,2,(np.argmax(p_save[2,2,:]['tid'],0))]['pos'])
# print((np.argmax(p_save[2,2,:]['tid'],0)))
# p_max[2,2]['tid']=(np.argmax(p_save[2,2,:]['tid'],0))
# print(p_max[2,2]['tid'])
# print(p_max[2,2])

for i in range(N):
    for j in range(3):
        p_max[i, j] = p_save[i, j, (np.argmax(p_save[i, j, :]['tid'], 0))]
        p_min[i, j] = p_save[i, j, (np.argmin(p_save[i, j, :]['tid'], 0))]

# for i in range(N):
    # print(i)
    # for j in range(3):
        # print(j)
        # print("max "+str(p_max[i, j]['pos']))
        # print("min "+str(p_min[i, j]['pos']))

p_x = p_save[:, 0, :]
p_y = p_save[:, 1, :]
p_z = p_save[:, 2, :]
p_xy = np.stack((p_x, p_y), 1)
# p_yx = np.stack((p_y,p_x),1)
# p_xz = np.stack((p_x,p_z),1)
# p_yz = np.stack((p_y,p_z),1)
# p_zx = np.stack((p_z,p_x),1)
# p_zy = np.stack((p_z,p_y),1)
rita.setup(p_xy, N, Nt, KE_save, PE_save, dt)
# rita.setup(p_xz, N, Nt, KE_save, PE_save, dt)

# v.40 fel med hastigheten i Y led för kropp 1, kanske är acceleratioenen som inte ändras? fixat, tror jag.
# v.41 ändrat så att allt blir float och kontrollerat att förhållandet mellan position och den resulterande
# accelerationen och positionen funkar på sammasätt i alla tre dimensioner.
# v.42 Ändrade datatyp från python float till numpy.longfloat som har ett större antal decimaler den kan hantera, detta
# fixade att förstora värden inte ändrades när programet kördes. Gjorde så att den bara skriver ut positons och
# hastighets värden ibland.
# v.43 måndag, Undersöka vilka variabler som måste vara longfloat för att kunna optimera programmet. Massan verkar
# påverka allting på ett rimligt sätt även om talen är har en 2-3 gånger så stor exponent. Fixat så att det går att
# slumpa fram start värden har fixat så att det är rätt sign på alla värden, alltså positiva massor och positiva eller
# negativa hastigheter eller positioner i alla tre dimensioner.
# tisdag, Säkerställt att alla datatyper är rätt, implementering av en enkel meny för att välja typ av simulation och
# hur länge den sak köra. Fixat kod för att få fram kinetisk och potentiel energi. Kontrollerat att energi
# uträkningarna blir rätt, men de blir fel om avståndet mellan kroppar blir för liten. Upptäckte att main loppen
# hämtade accelerationen flera gånger varje tidsteg då jag skulle printa energin, fixat nu. Gjort så att angivna
# tidvärden är mer intuativa nu.
# Onsdag, framställt bättre start värden för solen, jorden och månen utifrån jordens omlopps plan kring solen. Testat
# acceleration uträkningarna. Gjort så att positonerna och energin sparas.
# Torsdag, fixat så att det räknar ut max och min postitioner och skillnaden med vid vilket tidsteg sparat. Förenklat
# uträkningen av max och min position och tid väldigt mycket genom att ta bort onödiga arrayer. Skapat ett system för
# att få ut positionerna i bara två koordinater. Fixat datatyp på de sparade datatyper.
# Fredag, gjort så att ritningen av planeter och omloppsbanor använder två av dimensionerna från simulations koden,
# arbetet med skalan på det som ritas, lagt till merkurius, venus och mars med så bra startvärden som vi kunde hitta,
# gjort så att information om tiden och energin i systemet visas i hörnet av ritningen av omloppsbanorna.
