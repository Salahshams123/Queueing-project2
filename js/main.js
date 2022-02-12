

function showDete()
{
    document.getElementById("det").style.display = 'block';
    document.getElementById("sto").style.display = 'none';
}
function showStoc()
{
    document.getElementById("sto").style.display = 'block';
    document.getElementById("det").style.display = 'none';
}
function getInput(input = "")
{
    if (input[1] === "/")
        return parseFloat(input[0] / parseFloat(input[2]));
        
    else
        return parseFloat(input);
}

//Determanistic
let calc = function (lamda, mu, k, t, m)
{
    let n, n2 = -1;
    let res = [];
    let AT = 1 / lamda; //arrival time funtion in the arrival rate
    let ST = 1 / mu;    // service time funtional in the service rate
    let ti;
    if (lamda > mu)
    {
        
        if (k == 0 || isNaN(k))
        {
            if (t < AT)
                n = 0;
            else
                n = Math.floor(lamda * t) - Math.floor((mu * t) - (mu / lamda));
            
        }
        else
        {
            
            // find the time of occurrence of the first balk (ti)
            for (let i = 1; i <= (k - (mu / lamda)) / (lamda - mu); i++)
            {
                
                let s = Math.floor(lamda * i);
                let f = Math.floor((mu * i) - (mu / lamda));
                let U = s - f;
                if (k == U)
                {
                    ti = i;
                    break;
                }
            }

            if (t < AT)
                n = 0;
            else if (t >= AT && t < ti)
                n = Math.floor(lamda * t) - Math.floor((mu * t) - (mu / lamda));
            else if (t >= ti)
            {
                n = k - 1;
                n2 = k - 2;
            }
        }
    }
    else
    {
        if (lamda == mu)
        {
            n = m;
        }
        else
        {
            
            ti = Math.floor((m) / (mu - lamda));
            for (let i = 0; i <= ti; i++)
            {
                let s = Math.floor(lamda * i);
                let f = Math.floor(mu * i);
                let U = m + s - f;

                if (U == 0)
                {
                    ti = i;
                    break;
                }
            }
            let Q = Math.trunc(lamda * ti);
            if ( t < ti)
                n = m + Math.floor(t * lamda) - Math.floor(t * mu);
                
            else if (t >= ti)
            {
                n = 0;
                n2 = 1;
            }
            
        }
    }
    res[0] = n;
    res[1] = n2;
    res[2] = ti;
    return res;
}

console.log(calc(1 / 3, 1, 5, 0, 7));
function calcD_n()
{
    let lamda = getInput(document.getElementById("lamda").value);
    let mu = getInput(document.getElementById("mu").value);
    let k = getInput(document.getElementById("k").value);
    let t = getInput(document.getElementById("t").value);
    let m = getInput(document.getElementById("m").value);
    
    let info = calc(lamda, mu, k, t, m);
        
        if (info[1] != -1)
            document.getElementById("n").textContent = `n(t) = ${info[0]} or ${info[1]}`;
        else
            document.getElementById("n").textContent = `n(t) = ${info[0]}`;
    
}
function calc_wq(){
    let n = getInput(document.getElementById("n_").value);
    let lamda = getInput(document.getElementById("lamda").value);
    let mu = getInput(document.getElementById("mu").value);
    let k = getInput(document.getElementById("k").value);
    let t = getInput(document.getElementById("t").value);
    let m = getInput(document.getElementById("m").value);
    let AT = 1 / lamda; //arrival time funtion in the arrival rate
    let ST = 1 / mu;    // service time funtional in the service rate
    let wq, wq2 = -1;
    let info = calc(lamda, mu, k, t, m);
    let ti = info[2];
    if (mu < lamda)
    {
        if (n == 0)
            wq = 0;
        else if (n < Math.trunc(lamda * ti))
            wq = (ST - AT) * (n - 1);
        else
        {
            wq = (ST - AT) * (lamda * ti - 2);
            wq2 = (ST - AT) * (lamda * ti - 3);
        }
    }
    else
    {
        if (lamda == mu)
            wq = (m - 1) * ST;
        else
        {
            if (n == 0)
                wq = (m - 1) / (2 * mu);
            else if (n < Math.trunc(lamda * ti))
                wq = (m - 1 + n) * ST - n * (AT);
            else if (n >= Math.trunc(lamda * ti))
                wq = 0;
        }
    }
    console.log("sss")
    if (wq2 != -1)
        document.getElementById("w").textContent = `Wq(n) = ${wq} or ${wq2}`;
    else
        document.getElementById("w").textContent = `Wq(n) = ${wq}`;
}
function draw()
{
    
    let labels_ = [];
    let data_ = [];
    let  ti;
    let lamda = getInput(document.getElementById("lamda").value);
    let mu = getInput(document.getElementById("mu").value);
    let k = getInput(document.getElementById("k").value);
    let m = getInput(document.getElementById("m").value);
    let info = calc(lamda, mu, k, 0, m);
    if (info[2] == undefined)
    {
        ti = 30;
    }
    else
    {
        ti = info[2];
    }

    for (let i = 0; i < ti; i++)
        labels_[i] = i;
    
    for (let i = 0; i < ti; i++)
        data_[i] = calc(lamda, mu, k, i, m)[0];
    const ctx = document.getElementById("myChart");
    const myChart = new Chart(ctx, {
    type: 'line',
    data: {
        labels: labels_,
        datasets: [{
            label: 'n(t) = ',
            data: data_,
            borderColor: 'rgba(255, 99, 132, 1)',
            fill: false,
            stepped: true,
        }]
    },

        options: {
            responsive: true,
            interaction: {
                intersect: false,
                axis: 'x'

            },
            plugins: {
                title: {
                    display: true,
                    text: (ctx) => "number of customers in a D/D/1/k-1 system",
                }
            },
            scales: {
            y: {
                
                
                ticks: {
                    stepSize: 1,
                    color: 'black'
                }
            }
        },
        }
});

}
//Stochastic

function _M_M_1()
{
    //M/M/1
    document.getElementById("s_k").style.display = "none";
    document.getElementById("s_c").style.display = "none";
    document.getElementById("b_1").style.display = "block";
    document.getElementById("b_2").style.display = "none";
    document.getElementById("b_3").style.display = "none";
    document.getElementById("b_4").style.display = "none";
    
}
function calcS_1(){
    let lamda = getInput(document.getElementById("_lamda").value);
    let mu = getInput(document.getElementById("_mu").value);

    document.getElementById("_l").textContent =`L: ${(lamda / (mu - lamda))}` ;
    document.getElementById("_lq").textContent =`Lq: ${(lamda * lamda) / (mu * (mu - lamda))}` ;
    document.getElementById("_w").textContent =`W: ${(1 / (mu - lamda))}` ;
    document.getElementById("_wq").textContent =`Wq: ${(lamda) / (mu * (mu - lamda))}` ;
}

function _M_M_1_K()
{
    document.getElementById("s_k").style.display = "block";
    document.getElementById("s_c").style.display = "none";
    document.getElementById("b_2").style.display = "block";
    document.getElementById("b_1").style.display = "none";
    document.getElementById("b_3").style.display = "none";
    document.getElementById("b_4").style.display = "none";
}
function calcS_2()
{
    let lamda = getInput(document.getElementById("_lamda").value);
    let mu = getInput(document.getElementById("_mu").value);
    let k = getInput(document.getElementById("_k").value);
    let p0, pk;
    let ρ = (lamda / mu);
    let u = Math.pow(ρ, k + 1);
    if (ρ != 1)
        p0 = (1 - ρ) / (1 - u);
    else if (ρ == 1)
        p0 = 1 / (k + 1);
    
    let u1 = Math.pow(ρ, k);
    if (ρ == 1)
        pk = 1 / (k + 1);
    else
        pk = u1 * p0;
    
    let l = ((1 - ((k + 1) * u1) + (k * u)) / ((1 - ρ) * (1 - u))) * ρ;
    if (ρ == 1)
        document.getElementById("_l").textContent = k / 2;
    else
        document.getElementById("_l").textContent =`L: ${l}`;
    document.getElementById("_lq").textContent = `Lq: ${l - ((lamda* (1 - pk)) / mu)}`;
    document.getElementById("_w").textContent = `W: ${l / (lamda * (1 - pk))}`;
    let w = l / (lamda * (1 - pk));
    document.getElementById("_wq").textContent =`Wq: ${w - (1 / mu)}` ;

}

function _M_M_C()
{
    document.getElementById("s_c").style.display = "block";
    document.getElementById("s_k").style.display = "none";
    document.getElementById("b_3").style.display = "block";
    document.getElementById("b_1").style.display = "none";
    document.getElementById("b_2").style.display = "none";
    document.getElementById("b_4").style.display = "none";
}
function calcS_3()
{
    let lamda = getInput(document.getElementById("_lamda").value);
    let mu = getInput(document.getElementById("_mu").value);
    let c = getInput(document.getElementById("_c").value);
    let ρ, r, p0, pn;
    r = lamda / mu;
    ρ = r / c;
    function fact(a)
    {
        if (a == 0)
            return 1;
        else
            return a * fact(a - 1);
    }

    let sum = 0
    for (let i = 0; i <= c - 1; i++)
        sum += (Math.pow(r, i) / fact(i));
    
    if (ρ < 1)
    {
        p0 = 1 / (sum + ((Math.pow(r, c) * c) / (fact(c) * (c - r))));
    }
    else
        p0 = 1 / (sum + (((1 / fact(c)) * Math.pow(r, c))) * ((c * mu) / ((c * mu) - lamda)));
    
    let x = (c * mu) - lamda;
    let lq = ((Math.pow(r, c) * mu * lamda) / (fact(c - 1) * Math.pow(x, 2))) * p0;
    let l = lq + r;
    let wq = lq / lamda;
    let w = wq + (1 / mu);

    document.getElementById("_l").textContent = ` L: ${l}`;
    document.getElementById("_lq").textContent = `Lq: ${lq}`;
    document.getElementById("_w").textContent = `W: ${w}`;
    document.getElementById("_wq").textContent =`Wq: ${wq}` ;
}
function _M_M_C_K()
{
    document.getElementById("s_c").style.display = "block";
    document.getElementById("s_k").style.display = "block";
    document.getElementById("b_4").style.display = "block";
    document.getElementById("b_1").style.display = "none";
    document.getElementById("b_2").style.display = "none";
    document.getElementById("b_3").style.display = "none";
}
function calcS_4()
{
    let lamda = getInput(document.getElementById("_lamda").value);
    let mu = getInput(document.getElementById("_mu").value);
    let c = getInput(document.getElementById("_c").value);
    let k = getInput(document.getElementById("_k").value);
    let ρ, r, p0, pk;
    r = lamda / mu;
    ρ = r / c;
    function fact(a)
    {
        if (a == 0)
            return 1;
        else
            return a * fact(a - 1);
    }

    let sum = 0
    for (let i = 0; i <= c - 1; i++)
        sum += (Math.pow(r, i) / fact(i));
    
    
    if (ρ == 1)
    {
        let x = Math.pow(r, c) / fact(c);
        let y = k - c + 1;
        let z = sum + (x * y);
        p0 = 1 / z;
    }
    else
    {
        let x = Math.pow(r, c) / fact(c);
        let y = (1 - Math.pow(ρ, (k - c + 1))) / (1 - ρ);
        let z = sum + (x * y);
        p0 = 1 / z;
    }    
    
    if (k >= 0 && k < c)
        pk = (Math.pow(r, k) / fact(k)) * p0;
    else
        pk = (Math.pow(r, k) / (Math.pow(c, (k - c)) * fact(c))) * p0;

    
    let sum2 = 0;
    for (let i = 0; i <= c - 1; i++)
        sum2 += (c - i) * (Math.pow(r, i) / fact(i));
    
    
    let x = (ρ * Math.pow(r, c) * p0) / (fact(c) * Math.pow((1 - ρ), 2));
    let y = 1 - Math.pow(ρ, (k - c + 1)) - (1 - ρ) * (k - c + 1) * Math.pow(ρ, k - c);
    let lq = x * y;

    let l = lq + c - (sum2 * p0);

    let wq = lq / (lamda * (1 - pk));

    let w = l / (lamda * (1 - pk));
    document.getElementById("_l").textContent = ` L: ${l}`;
    document.getElementById("_lq").textContent = `Lq: ${lq}`;
    document.getElementById("_w").textContent = `W: ${w}`;
    document.getElementById("_wq").textContent =`Wq: ${wq}` ;
}
