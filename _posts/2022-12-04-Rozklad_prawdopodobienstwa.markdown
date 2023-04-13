---
layout: post
title:  "Rozkład prawdopodobieństwa"
date:   2022-12-04 14:34:47 +0100
---

# Rzut monetą

W jednorazowym rzucie monetą (próba Bernouliego) prawdopodobieństwo sukcesu wynosi $p=1/2$ np. wyrzucenie reszki. Jeśli eksperyment będzie polegać na trzykrotnym rzucie monetą (trzy próby Bernouliego) to otrzymamy próbę losową z rozkładu Bernoulliego (szczególny przypadek rozkładu dwumianowego) gdzie $1$ oznacza sukces natomiast $0$ porażkę.

Próba losowa z rozkładu dwumianowego o parametrach $n=3$ oraz $p=0,5$ może oznaczać wynik:

$$
(O,O,R)\quad\mathrm{lub}\quad (0,0,1)
$$

Dodajmy, że wszystkich możliwych wariantów jest osiem ponieważ $2^3=8$:

$$
\{(O,O,O),(R,O,O),(O,R,O),(R,R,O),\
         (O,O,R),(R,O,R),(O,R,R),(R,R,R)\}
$$

Na podstawie powyższego zbioru można wyznaczyć liczbę sukcesów/reszek:

$$
0, 1, 1, 2, 1, 2, 2, 3
$$

```julia
using Combinatorics, IterTools, Distributions, StatsBase, Plots, LinearAlgebra

# wszystkie możliwe warianty:
M = product([0,1],[0,1],[0,1]);
sol = vec(collect.(M))
8-element Vector{Vector{Int64}}:
 [0, 0, 0]
 [1, 0, 0]
 [0, 1, 0]
 [1, 1, 0]
 [0, 0, 1]
 [1, 0, 1]
 [0, 1, 1]
 [1, 1, 1]
# liczba sukcesów tj. wypadnięcie reszki dla każdego wariantu:
res = sum.(sol)
8-element Vector{Int64}:
 0
 1
 1
 2
 1
 2
 2
 3
sort(countmap(res))
OrderedCollections.OrderedDict{Int64, Int64} with 4 entries:
  0 => 1
  1 => 3
  2 => 3
  3 => 1
```
Dla większej liczby rzutów monetą wyznaczenie wszystkich możliwych wyników może być kłopotliwe ponieważ ich liczba rośnie wykładniczo.
Rozwiązaniem może być wykorzystanie generatora liczb losowych z rozkładu Bernulliego aby wygenerować dowolną liczbę prób Bernoulliego np. $10000$. Poniżej symulacja rozkładu liczby sukcesów w trzech rzutach monetą.
```julia
# 10000 prób losowych Bernoulliego:
S = rand(Binomial(1,0.5),3,10000)
3×10000 Matrix{Int64}:
 1  1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  …  1  1  1  0  0  0  0  0  1
 0  0  1  0  0  0  1  0  0  1  0  1  0  0  1  1     0  1  1  0  0  0  0  1  0
 0  0  0  0  0  1  0  1  1  1  1  1  1  1  1  0     1  0  1  0  0  0  0  0  0

res1 = vec(sum(S, dims=1));
sort(countmap(res1))
OrderedCollections.OrderedDict{Int64, Int64} with 4 entries:
  0 => 1264
  1 => 3726
  2 => 3756
  3 => 1254
```
Jeśli liczbą sukcesów będzie to, że reszka wystąpi dokładnie dwa razy (są trzy takie przypadki)
to dokładne prawdopodobieństwo tego zdarzenia wyniesie:

$$
P(X=2)=3/8=0.375
$$

```julia
# dokładne prawdopodobieństwo:
mean(res.==2)
0.375
# Binomial prawdopodobieństwo:
pdf.(Binomial(3,0.5),2)
0.3750000000000001
# symulacja dokładnego prawdopodobieństwa p-value=3756/10000:
mean(res1.==2)
0.3756
```

Do graficznej prezentacji dokładnego rozkładu prawdopodobieństwa liczby sukcesów np. trzykrotnego rzutu monetą można wykorzystać rozkład dwumianowy $\mathrm{B}(3,1/2)$.
```julia
X, Y = 0:3, round.(pdf.(Binomial(3,0.5),0:3),digits=3);
bar(X,Y,labels="Binomial",
    fillcolor="blue",fillalpha=0.3,lc="white",
    xlab="Liczba sukcesów", ylab="Prawdopodobieństwo",
    titlefontsize=10, title="Dwumianowy rozkład prawdopodobieństwa")
annotate!(X,Y,Y,8)
savefig("/home/krz/Pulpit/blog/figure/a01.png")
```

![a01](/assets/a01.png)

![a01]({{ "/assets/a01.png" | relative_url }}) 

<figure>
<center>
<img alt="png" src="/figure/a01.png">
</center>
</figure>

Prawdopodobieństwo tego, że np. reszka wystąpi co najwyżej dwa razy (skumulowane prawdopodobieństwo) będzie równe:

$$P(X\leq 2)=P(X=2)+P(X=1)+P(X=0)=0.875$$

```julia
# dokładne prawdopodobieństwo:
mean(res.<=2)
0.875
# Binomial prawdopodobieństwo:
cdf.(Binomial(3,0.5),2)
0.875
# symulacja dokładnego prawdopodobieństwa:
mean(res1.<=2)
0.8746
```

Skumulowane prawdopodobieństwo:

```julia
X, Y = 0:3, round.(cdf.(Binomial(3,0.5),0:3),digits=3);
bar(X,Y,labels="Binomial",legend=:topleft,
    fillcolor="blue",fillalpha=0.3,lc="white",
    xlab="Liczba sukcesów", ylab="Skumulowane prawdopodobieństwo",
    titlefontsize=10, title="Skumulowany rozkład prawdopodobieństwa")
annotate!(X,Y,Y,8)
savefig("/home/krz/Pulpit/blog/figure/a02.png")
```
<figure>
<center>
<img alt="png" src="/figure/a02.png">
</center>
</figure>

# Rzut kostką do gry

W jednorazowym rzucie kostką do gry prawdopodobieństwo sukcesu wyrzucenia określonej liczby oczek wynosi $p=1/6$.
Warto zauważyć, że liczba wszystkich możliwych wariantów rośnie wykładniczo $6^n$. Zatem dla dwóch rzutów kostką do gry wszystkich możliwoścu będzie $6^2$.

Na przykładzie dwukrotnego rzutu kostką do gry można wyznaczyć wszystkie możliwe wyrzucone pary oczek oraz zaznaczyć
warianty w których suma oczek jest równa $6$.

$$
\begin{Bmatrix}
    \textit{(1,1)} & \textit{(1,2)} & \textit{(1,3)} & \textit{(1,4)} & \textbf{(1,5)} & \textit{(1,6)}\\
    \textit{(2,1)} & \textit{(2,2)} & \textit{(2,3)} & \textbf{(2,4)} & \textit{(2,5)} & \textit{(2,6)}\\
    \textit{(3,1)} & \textit{(3,2)} & \textbf{(3,3)} & \textit{(3,4)} & \textit{(3,5)} & \textit{(3,6)}\\
    \textit{(4,1)} & \textbf{(4,2)} & \textit{(4,3)} & \textit{(4,4)} & \textit{(4,5)} & \textit{(4,6)}\\
    \textbf{(5,1)} & \textit{(5,2)} & \textit{(5,3)} & \textit{(5,4)} & \textit{(5,5)} & \textit{(5,6)}\\
    \textit{(6,1)} & \textit{(6,2)} & \textit{(6,3)} & \textit{(6,4)} & \textit{(6,5)} & \textit{(6,6)}
\end{Bmatrix}
$$

Ponieważ ilość wszystkich możliwości jest równa 36 to prawdopodobieństwo tego, że
suma oczek będzie równa 6 wynosi:

$$P(X=6)=5/36=0.1389$$

Dla dużych $n$ wyznaczenie dokładnego rozkładu prawdopodobieństwa może być utrudnione ze względu na dużą ilość obliczeń. Rozwiązaniem może być aproksymacja czyli przybliżanie rozkładu dokładnego za pomocą symulacji. 
```julia
# symulacja dla ośmiu rzutów kostką:
m = [sum(sample(1:6,8)) for i in 1:100000];
# parametry rozkładu normalnego:
mu1, sd1 = mean(m), std(m)
(28.01644, 4.820769868456823)
```

Rozkład normalny $\mathrm{N(\mu,\sigma)}$ w którym:

$$\mu=np,\quad \sigma=\sqrt{np(1 − p)}$$

Liczba sukcesów $n_1$ ma rozkład dwumianowy $\mathrm{Bin}(n,p)$ o parametrach:

$$n=-\mu^2/(\sigma^2-\mu),\quad p=-(\sigma^2-\mu)/\mu$$

```julia
# parametry rozkładu dwumianowego:
n1, p1 = (-mu1^2/(sd1^2-mu1), -(sd1^2-mu1)/mu1)
(164.32566530379117, 0.17049339157219057)
```

{::options parse_block_html="true" /}
<details><summary markdown="span">**CODE**</summary>
```julia
c1 = countmap(m);
bar(collect(keys(c1)),values(c1)./length(m),labels="random Binomial",
    fillcolor="blue",fillalpha=0.3,lc="white",
    titlefontsize=10, title="Rozkład prawdopodobieństwa liczby sukcesów")
plot!(x -> pdf(Normal(mu1,sd1),x), color=:green, alpha=0.45, lw=3,
      labels="Normal")
scatter!(collect(keys(c1)),pdf.(Binomial(round(n1), p1),collect(keys(c1))),labels="Binomial",
         markerstrokecolor = :red, markercolor = :red, alpha = 0.5)
savefig("/home/krz/Pulpit/blog/figure/a03.png")
```
</details>
<br/>
{::options parse_block_html="false" /}

<figure>
<center>
<img alt="png" src="/figure/a03.png">
</center>
</figure>

Dla dużych $n$ prawdopodobieństwo sukcesu $p=n_1/n$ ma asymptotycznie rozkład normalny gdzie:

$$\mu=p,\quad \sigma=\sqrt{p(1-p)/n}$$

Dodatkowo można założyć, że proporcja sukcesu $p$ ma rozkład beta $\mathrm{Bet}(n_1,n_0)$ o parametrach $n_1$ oraz $n_0$ gdzie:

$$n_1=\mu\cdot\eta,\quad n_0=(1-\mu)\cdot\eta$$

dla $\eta=\frac{(1-\mu)\mu}{\sigma^2}-1$.

```julia
# symulacja:
M = m./48;
# parametry rozkładu normalnego:
MU1, SD1 = mean(M), std(M)
(0.5836758333333333, 0.10043270559285054)
# parametry rozkładu beta:
N = ((MU1*(1-MU1))/SD1^2)-1;
a, b  = MU1*N, (1-MU1)*N
(13.47759977009273, 9.613299322170636)
```

{::options parse_block_html="true" /}
<details><summary markdown="span">**CODE**</summary>
```julia
plot(normalize(fit(Histogram, M, closed=:left, nbins=30), mode=:pdf),
     fillcolor="blue",fillalpha=0.3,lc="white",labels="random",
     titlefontsize=10, title="Rozkład prawdopodobieństwa proporcji sukcesów")
plot!(x -> pdf(Normal(MU1,SD1),x), color=:green, alpha=0.45, lw=3,
      labels="Normal")
plot!(x -> pdf(Beta(a,b),x), color=:red, alpha=0.4, lw=3,
      labels="Beta")
savefig("/home/krz/Pulpit/blog/figure/a04.png")
```
</details>
<br/>
{::options parse_block_html="false" /}

<figure>
<center>
<img alt="png" src="/figure/a04.png">
</center>
</figure>

Za pomocą rozkładu normalnego można aproksymować (przybliżać) rozkłady dyskretne oraz ciągłe. W przypadku przybliżeń rozkładów dyskretnych należy pamiętać o stosowaniu korekty na ciągłość:

$$
P(x_1<X<x_2)=\int_{x_1-0,5}^{x_2+0,5}f(x\,|\,\theta)\,dx
$$

Poniżej przykład wyznaczenia prawdopodobieństwa $P(X\leq 15)$:

```julia
# empiryczna dystrybuanta:
ecdf(m)(15+0.5)
0.00366
# dystrybuanta rozkładu normalnego:
cdf(Normal(mu1,sd1),15+0.5)
0.004710904600789341
# dystrybuanta rozkładu dwumianowego:
cdf(Binomial(round(n1), p1),15)
0.00292650828913315
```

