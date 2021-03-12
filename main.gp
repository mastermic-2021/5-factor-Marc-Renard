encode(m)=fromdigits(Vec(Vecsmall(m)),128);
decode(m)=Strchr(digits(m, 128));

[n,c] = readvec("input.txt");

dechiffrer(c,p,q)={   \\déchiffrement d'un chiffré par le système de Rabin
	chiffre=Mod(c[1],n); \\récupération de m²mod(n)
	kornecker=c[2];
	mod=c[3];
	modp=chiffre^((p+1)/4); \\calcul de la racine carré mod p
	modq=chiffre^((q+1)/4);  \\calcul de la racine carré mod q
	[u,v]=bezout(p,q);
	\\on a pu+vq=1 donc pu = 1 mod q et donc pu*modq = modq mod(q)
	\\Même chose avec pv
	m1=(u*p*modq+v*q*modp);
	m2=n-m1;
	m3=(u*p*modq-v*q*modp);
	m4=n-m3;
	mpot=[lift(m1),lift(m2),lift(m3),lift(m4)];
	for(i=1,4,if((kronecker(mpot[i],n)==c[2])&&(mpot[i]%2==c[3]),return(mpot[i])));
}

PMUPollard(n) = {
	i=2;
	a=Mod(2,n);
	\\calcul de a^(i!) jusqu'à ce qui i! soit un multiple de p-1 et donc que pgcd(a^(i!)-1,n)>1
	while(1,a=a^i;d=gcd(lift(a-1),n);if(d>1,break);i++);
	d;
};

p=PMUPollard(n);
a=dechiffrer(c,p,n/p);
a=lift(a);
message=decode(a);
print(message);

