var trig=0; //trigonometric function? or exponential? used to restrict domain of calculations when optimization is needed

class Function {
  //var Roots=[];
  constructor(degree, a) {
    this.co = [];
    this.deg = degree;
    this.roots = [];
    this.rootsNum=0;
    this.forbiddenArg = []; //shall remain empty
    this.visibleRoots=[];
    this.domainFound=0;
    for (let i = 0; i <= this.deg; i++) {
      if (typeof a[i] !== 'undefined')
        this.co[i] = a[i];
      else
        this.co[i] = 0;
    }
  }
  Val(x) {
    if(this.isInDomain(x)==1){
        let s = 0;
        let x_pow = 1;
        for (let i = this.deg; i >= 0; i--) {
          if (typeof this.co[i] === 'undefined')
            this.co[i] = 0;
          s += (this.co[i] * x_pow);
          if (i != 0)
            x_pow = x_pow * x;
        }
        return s;
    }
    else {
        return 0;
    }      
  }
    
  Write() {
    let text = "";
    for (let i = 0; i <= this.deg; i++) {
      if (this.co[i] != 0) {
        if (i == 0) {
          if (this.co[i] != 1 && this.co[i] != -1)
            text += this.co[i];
          else if (this.co[i] == -1)
            text += '-';
        }
        else if (this.co[i] > 0) {
          if (this.co[i] != 1 || i == this.deg)
            text += '+' + this.co[i];
          else
            text += '+';
        }
        else if (this.co[i] != -1)
          text += this.co[i];
        else
          text += "-";
        if (i < this.deg - 1)
          text += "x^" + (this.deg - i);
        else if (i < this.deg)
          text += "x";
      }
    }
    return text;
  }
  Integrate(a, b) {
    if (this.deg == 0)
      return (this.Val(a) * (b - a));
    else if (this.deg == 1)
      return ((this.Val(a) + this.Val(b)) * (b - a) / 2);
    else {
      let n = 50;
      let s = 0;
      let h = (b - a) / (4 * n);
      for (let i = 0; i < 2 * n; i++) {
        s += h / 3 * (this.Val(a + 2 * i * h) + this.Val(a + 2 * (i + 1) * h) + 4 * this.Val(a + h * (1 + 2 * i)));
      }
      return s;
    }
  }
  Differentiate(a) {
    let h = 0.000000000001;
    return (this.Val(a + h) - this.Val(a)) / h;
  }
  DifferentiateN(a, n) {
    let h;
    if (n == 1)
      return this.Differentiate(a);
    else {
      h = 0.00001 * this.DifferentiateN(a, n - 1);
      if (h == 0)
        h = 0.000000000001
      return (this.DifferentiateN(a + h, n - 1) - this.DifferentiateN(a, n - 1)) / h;
    }
  }
  DrawTangent(x) {
    let a = this.Differentiate(x);
    let b = this.Val(x) - a * x;
    let G = new Function(1, [a, b]);
    return G;
  }
  FindRoots() {
    //ITP method
    let R = []; //roots
    let d = 0; //number of roots
    let Ri = 0; //root iterator
    //declaration
    let eps = 0.001 //precision
    let a = -1000000; //we start in (relative) minus infinity
    if(this.isInDomain(a)==0)
    {
        a = this.moveToDomain(a);
    }
    let b = a; //upper boundary, iterated
    let k1 = 10; //scaling truncation (0;inf)
    let k2 = 2; //exponential part of truncation [1;1+phi)
    let n0 = 10; //projection scalar, unsure of role [0;inf)
    let xm; //x1/2
    let xf; //regula falsi
    let nm; //log
    let x; //sought x
    //help variables
    let fa = this.Val(a); //f(a)
    let fb; //f(b)
    let fx; //f(x)
    let sig; //sign in truncation
    let trunc; //estimator vector base in truncation
    let xt; //truncated x
    let proj; //estimator in projection
    let j = 0; //iterator
    let flag = 0; //any roots?
    let lR = 0; //last root calculated
    //b iteration
    while (true) {
      while (b <= 1000000 && Math.sign(fa) == Math.sign(this.Val(b))) {
        a = b;
        fa = this.Val(a);
        b += 0.2;
        if(this.isInDomain(b)==0)
            b = this.moveToDomain(b);
      }
      if (b > 1000000 && flag == 0)
        return "No real roots";
      else if (b > 100000) {
        this.roots = R;
        this.rootsNum = Ri;
        return R;
      }
      flag = 1;
      j = 0;
      if (this.Val(b) == 0) {
        a = b;
        x = b;
      }
      while (Math.abs(b - a) > 2 * eps) {
        fb = this.Val(b);
        nm = Math.ceil(Math.log2((b - a) / (2 * eps)));
        //interpolation
        xm = (a + b) / 2;
        xf = (b * fa - a * fb) / (fa - fb);
        //truncation
        sig = Math.sign(xm - xf);
        trunc = Math.min(sig * (xm - xf), k1 * Math.pow(Math.abs(b - a), k2));
        xt = xf + sig * trunc;
        //projection
        proj = Math.min(eps * Math.pow(2, nm + n0 - j) - (b - a) / 2, Math.abs(xt - xm));
        x = xm - sig * proj;
        fx = this.Val(x);
        if (fx == 0) {
          a = x;
          b = x;
        }
        else if (Math.sign(fx) == Math.sign(fa))
          a = x;
        else
          b = x;
        fa = this.Val(a);
        j++;
        if (j > 50)
          return 'undefined';
      }
      if (Ri == 0 || Math.abs((a + b) / 2 - lR) > 10 * eps) {
        R[Ri] = (a + b) / 2;
        lR = R[Ri];
        Ri++;
      }
      a = R[Ri - 1] + 10 * eps;
      b = a + 10 * eps;
      fa = this.Val(a);
      fb = this.Val(b);
    }

  }
  Merge(G) {
    function Reindex(F, dg, dl)
    {
      let b = [];
      for (let i = 0; i < (dg - dl); i++)
        b[i] = 0;
      for (let i = dg - dl; i <= dg; i++)
        b[i] = F.co[i - dg + dl];
      return b;
    }
    if (this.deg > G.deg)
      G.co = Reindex(G, this.deg, G.deg);
    else if (this.deg < G.deg)
      this.co = Reindex(this, G.deg, this.deg);
    this.deg = Math.max(this.deg, G.deg);
    for (let i = 0; i <= this.deg; i++) {
      if (typeof this.co[i] === 'undefined')
        this.co[i] = 0;
      if (typeof G.co[i] === 'undefined')
        G.co[i] = 0;
      this.co[i] += G.co[i];
    }
  }
  FindMaximum(a, b) {
    if(this.isInDomain(a)==0)
        a=this.moveToDomain(a);
    let dx = Math.min(0.1, (b - a) / 50);
    let fmax = this.Val(a);
    let fcur;
    for (let x = a; x <= b; x += dx) {
      if(this.isInDomain(x)==1){
          fcur = this.Val(x);
          if (fcur > fmax)
            fmax = fcur;
      }
    }
    if(this.isInDomain(b)==1){
        if (this.Val(b) > fmax)
            fmax = this.Val(b);
    }
    return fmax;
  }
  FindMinimum(a, b) {
    if(this.isInDomain(a)==0)
        a=this.moveToDomain(a);
    let dx = Math.min(0.1, (b - a) / 50);
    let fmin = this.Val(a);
    let fcur;
    for (let x = a; x <= b; x += dx) {
      if(this.isInDomain(x)==1){
        fcur = this.Val(x);
        if (fcur < fmin)
        fmin = fcur;
      }
    }
    if(this.isInDomain(b)==1){
        if (this.Val(b) < fmin)
            fmin = this.Val(b);
    }
    return fmin;
  }
    
  findDomain()
  {
      this.domainFound=1;
  }
    
  isInDomain(x)
  {
      let is=1;
      this.forbiddenArg.forEach(element => {if(x>=element[0]&&x<=element[1]) {is=0;}})
      return is;
  }
  moveToDomain(x)
  {
      let lastForbidden;
      this.forbiddenArg.forEach(element => {if(x>=element[0]&&x<=element[1]) lastForbidden=element[1];})
      let soughtX=lastForbidden+0.1;
      if(this.isInDomain(soughtX)==0)
          this.moveToDomain(soughtX);
      return soughtX;
  }
  remapDomain(left,right){
      let firstGood;
      let lastGood;
      let nonEmpty=0;
      let i=0;
      let tempList=[];
      if(this.forbiddenArg.length>0 && typeof this.forbiddenArg[0]!=='undefined'){
          while(this.forbiddenArg[i][1]<left&& typeof this.forbiddenArg[i+1]!== 'undefined' ){
              i++;
          }
          if(this.forbiddenArg[i][1]>=left){
              firstGood=i;
              nonEmpty=1;
          }
          while(this.forbiddenArg[i][1]<right&& typeof this.forbiddenArg[i+1]!== 'undefined' ){
              i++;
          }
          if(this.forbiddenArg[i][0]<=right){
              lastGood=i;
          }
          else if(i-1!=firstGood)
              lastGood=i-1;
          else
              nonEmpty=0;
      }
      if(nonEmpty==1){
          for(i=firstGood;i<=lastGood;i++)
              tempList.push(this.forbiddenArg[i]);
          this.forbiddenArg=[];
          this.forbiddenArg.push([-1000000,left-0.5]);
          tempList.forEach(element => this.forbiddenArg.push(element));
          this.forbiddenArg.push([right+0.5,1000000]);
      }
      else{
          this.forbiddenArg=[];
          this.forbiddenArg.push([-1000000,left-0.5]);
          this.forbiddenArg.push([right+0.5,1000000]);
      }
      
  }

}

class CompositeFunction extends Function {
    
    constructor(functionA,functionB,compositionType)
    {
        super();
        this.funcType=compositionType;
        this.funcA=functionA;
        this.funcB=functionB;
        this.forbiddenArg = [];
        this.domainFound = 0;
        this.rootsNum;
        if(this.funcType=="sin" || this.funcType=="cos" || this.funcType=="tan" || this.funcType=="pow"){
            trig=1;
        }
    }
    
    findDomain(){
        this.domainFound=1;
        switch(this.funcType){
            case "frac":
                this.funcB.FindRoots();
                this.funcB.roots.forEach(element => this.configureForbiddenArg(element,0.075));
                break;
            //case "sin": nothing
            //case "cos": nothing
            case "tan":
                let tempFunc = new CompositeFunction(this.funcA,this.funcB,"cos");
                tempFunc.FindRoots();
                tempFunc.roots.forEach(element => {this.configureForbiddenArg(element,0.1)});
                break;
            case "log":
                this.findNegative(this.funcB);
                break;
            case "ln":
                this.findNegative(this.funcB);
                break;
            case "pow": //it omits the problem of 0^0 as for now (I think)
                if(this.funcA.deg>0||!Number.isInteger(this.funcA.co[0]))
                    this.findNegative(this.funcB);
                break;
            case "sqrt":
                if(!Number.isInteger(this.funcA.co[0])||this.funcA.co[0]%2==0)
                    this.findNegative(this.funcB);
                break;
            //case "abs": nothing
            case "arcsin":
                let temp1 = new Function(0,[1]);
                let temp2= new Function(0,[-1]);
                let tempFunc1 = new HolisticFunction([temp1,this.funcB]);
                let tempFunc2 = new HolisticFunction([temp2,this.funcB]);
                this.findNegative(tempFunc1); //[-1;inf)
                this.findPositive(tempFunc2); //(-inf;1]
                break;
            case "arccos":
                let tempb1 = new Function(0,[1]);
                let tempb2= new Function(0,[-1]);
                let tempFuncb1 = new HolisticFunction([tempb1,this.funcB]);
                let tempFuncb2 = new HolisticFunction([tempb2,this.funcB]);
                this.findNegative(tempFuncb1); //[-1;inf)
                this.findPositive(tempFuncb2); //(-inf;1]
                break;
            case "product":
                if(this.funcA.domainFound==0)
                    this.funcA.findDomain();
                if(this.funcB.domainFound==0)
                    this.funcB.findDomain();
                this.forbiddenArg=this.funcA.forbiddenArg.concat(this.funcB.forbiddenArg);
                break;
            //case "arctan": nothing
                
                
        }
    }
    
    findNegative(func)
    {
        func.FindRoots();
        let x=-1000000;
        let sig=Math.sign(func.Val(x));
        let i=0;
        while(i<func.rootsNum){
            if(sig==-1)
            {
                this.configureForbiddenArgRange(x,func.roots[i],0.05);
                i+=2;
            }
            else
                i++;
            x=func.roots[i-1];
            sig=-1;
        }
        if(i==func.rootsNum)
            if(sig==-1)
                this.configureForbiddenArgRange(x,1000000,0.05);
    }
    findPositive(func)
    {
        func.FindRoots();
        let x=-1000000;
        let sig=Math.sign(func.Val(x));
        let i=0;
        while(i<func.rootsNum){
            if(sig==1)
            {
                this.configureForbiddenArgRange(x,func.roots[i],0.05);
                i+=2;
            }
            else
                i++;
            x=func.roots[i-1];
            sig=1;
        }
        if(i==func.rootsNum)
            if(sig==1)
                this.configureForbiddenArgRange(x,1000000,0.05);
    }
    
    configureForbiddenArg(forbidden,extension)
    {
        this.forbiddenArg.push([forbidden-extension,forbidden+extension]);
    }
    configureForbiddenArgRange(forbidden1,forbidden2,extension)
    {
        this.forbiddenArg.push([forbidden1-extension,forbidden2+extension]);
    }
    
    Val(x) {
        if(!this.domainFound)
            this.findDomain();
        switch(this.funcType){
            case "frac":
                return (this.funcA.Val(x)/this.funcB.Val(x));
                break;
            case "sin": //for trigonometry funcA defines factor a in a*function(x)
                return this.funcA.co[0]*Math.sin(this.funcB.Val(x));
            case "cos":
                return this.funcA.co[0]*Math.cos(this.funcB.Val(x));
            case "tan":
                return this.funcA.co[0]*Math.tan(this.funcB.Val(x));
            case "log": //for logarithms funcA defines factor a in a*function(x), then the base of the logarithm
                return this.funcA.co[1]*Math.log(this.funcB.Val(x))/Math.log(this.funcA.co[0]);
            case "ln": //for logarithms funcA defines factor a in a*function(x)
                return this.funcA.co[0]*Math.log(this.funcB.Val(x));
            case "pow":
                return Math.pow(this.funcB.Val(x),this.funcA.Val(x));
            case "sqrt": //funcA defines factor a in a*function(x), then the degree of the root
                return this.funcA.co[1]*Math.pow(this.funcB.Val(x),1/this.funcA.co[0])
            case "abs":
                return this.funcA.co[0]*Math.abs(this.funcB.Val(x));
            case "arcsin":
                return this.funcA.co[0]*Math.asin(this.funcB.Val(x));
            case "arccos":
                return this.funcA.co[0]*Math.acos(this.funcB.Val(x));
            case "arctan":
                return this.funcA.co[0]*Math.atan(this.funcB.Val(x));
            case "product":
                return this.funcA.Val(x)*this.funcB.Val(x);
        }
    }
    
    Write() {
        let text="";
        switch(this.funcType){
                case "frac":
                    text+="\\frac{"+this.funcA.Write()+"}{"+this.funcB.Write()+"}";
                    break;
        }
        return text;
    }
    
}

class HolisticFunction extends CompositeFunction {
    constructor(functions){
        super();
        this.funcArray=[];
        if(typeof functions!==undefined)
            this.funcArray=functions;
        this.forbiddenArg=[];
        this.deg=0; //might be needed for some functions to execute properly, perhaps (e.g. integrate)
        this.domainFound=0;
        this.rootsNum;
    }

    Val(x){
        if(!this.domainFound)
        {
            this.findDomain();
        }
        let s=0;
        this.funcArray.forEach(element => {s+=element.Val(x)})
        return s;
    }
    
    Merge(G){
        if(G.constructor.name!="HolisticFunction")
            this.funcArray.push(G);
        else
            this.funcArray.push(G.funcArray); //perhaps unneeded, shall see
    }
    
    findDomain(){
        this.domainFound=1;
        this.funcArray.forEach(element => {element.findDomain(); element.forbiddenArg.forEach(subelement => this.forbiddenArg.push(subelement)); })
        this.funcArray.forEach(element => {element.forbiddenArg=this.forbiddenArg;})
    }
    
}

function Process(eq, x, y, t,mode) {
  let i = 0;
  let a; //coefficient
  let p; //power
  let k = 0;
  let flag; //0th power?
  let polynomialFlag=0; //Function class needed to merge?
  let bracket=0; //[]?
  let factor = 1; //1 or -1
  let funcType; //CompositeFunction type
  let f1temp;
  //let f2temp; //perhaps unneeded, shall see
  let hol1 = [];
  let hol2 = [];
  //fArray = []; //Function array
  //compArray = []; //CompositeFunction array
  let F = new Function();
  let Func = new HolisticFunction([]); //our main function
  if(mode=="ultimate")
      trig=0;
  while (eq.length > 0) {
    i = 0;
    flag = 0; //0th power?
    factor = 1;
    if (eq[0] == '-') {
      factor = -1;
      eq = eq.slice(1, eq.length);
    }
    while (typeof eq[i] !== 'undefined' && eq[i] != 'x' && eq[i]!='\\' && !((eq[i] == '+' || eq[i] == '-') && i != 0)) //iteruje się aż do znalezienia końca współczynnika
    {
      i++;
    }
    if (eq[i] != 'x')
      flag = 1;
    if (i == 1 && eq[0] == '+')
      a = 1;
    else if (i == 0)
      a = 1 * factor;
    else
      a = factor * Number(eq.slice(0, i));
      
    
    if(eq[i]=='\\'){
        Func.deg=2; //might be needed for integral to execute properly
        eq = eq.slice(i+1, eq.length); //now we cut \
        i=0;
        while (eq[i]!='[' && eq[i]!='{') //iteruje się aż do znalezienia końca formuły
        {
            i++;
        }
        funcType=eq.slice(0,i);
        if(eq[i]=='[')
            bracket=1;
        eq = eq.slice(i+1, eq.length); //now we cut formula
        i=0;
        if(bracket)
        {
            while(eq[i]!=']'){
                i++;
            }
            f1temp=Number(eq.slice(0, i));
            eq = eq.slice(i+2, eq.length); //now we cut to the content of the second bracket
            i=0;
        }
        while(eq[i]!='}'&&typeof eq[i]!=='undefined')
        {
            i++;
        }
        
        hol1.push(Process(eq.slice(0,i),0,0,0,"temp"));
        
        eq = eq.slice(i+1,eq.length);
        
        i=0;
        
        if(eq[0]=='}')
            eq=eq.slice(i+1,eq.length);
        
        if(eq[0]=='{')
        {
            eq = eq.slice(1,eq.length);
            while(eq[i]!='}'&&typeof eq[i]!=='undefined')
            {
                i++;
            }
            hol2.push(Process(eq.slice(0,i),0,0,0,"temp"));
            eq = eq.slice(i+1,eq.length);
            i=0;
            if(funcType!="pow") //to-do
                Func.funcArray.push(new CompositeFunction(hol1.slice(-1)[0],hol2.slice(-1)[0],funcType));
            else{
                Func.funcArray.push(new CompositeFunction(new Function(0,[a]),new CompositeFunction(hol1.slice(-1)[0],hol2.slice(-1)[0],"pow"),"product"));
            }
        }
        else
        {
            if(bracket)
                if(funcType!="pow")
                    Func.funcArray.push(new CompositeFunction(new Function(1,[f1temp,a]),hol1.slice(-1)[0],funcType));
                else
                    Func.funcArray.push(new CompositeFunction(new Function(0,[a]),new CompositeFunction(new Function(0,[f1temp]),hol1.slice(-1)[0],"pow"),"product"));
            else
                Func.funcArray.push(new CompositeFunction(new Function(0,[a]),hol1.slice(-1)[0],funcType));
        }
        bracket=0;
        
    }
    //composite function ends here
      
    else{
        polynomialFlag=1;

        eq = eq.slice(Math.min(i, eq.length - 1), eq.length);
        if (eq[0] == 'x')
          eq = eq.slice(1, eq.length);
        if (eq[0] == '^')
          eq = eq.slice(1, eq.length);
        i = 0;
        while (typeof eq[i] !== 'undefined' && eq[i] != 'x' && eq[i] != '+' && eq[i] != '-') {
          i++;
        }
        if (flag)
          p = 0;
        else if (i == 0)
          p = 1;
        else
          p = Number(eq.slice(0, i));
        if (k == 0) {
          F.deg = p;
          F.co[0] = a;
        }
        else {
          G = new Function(p, [0]);
          G.co[0] = a;
          F.Merge(G);
        }
        k = 1;
        eq = eq.slice(i, eq.length);
    }
  }
  if(polynomialFlag)
      Func.Merge(F);
  if(mode=="ultimate")
  {
      let A = Func.DrawTangent(t);
      if(trig==1)
          InitTrig(Func,x,y,t);
      document.getElementById("IntegrationAnswer").innerHTML = Func.Integrate(x, y);
      document.getElementById("1DerivativeAnswer").innerHTML = Func.Differentiate(t);
      //document.getElementById("2DerivativeAnswer").innerHTML = Func.DifferentiateN(t, 2);
      document.getElementById("TangentAnswer").innerHTML = A.Write();
      Func.FindRoots();
      DrawGraph(Func, A, x, y, t);
      if(Func.rootsNum>0)
        document.getElementById("RootsAnswer").innerHTML = Func.visibleRoots;
      else
        document.getElementById("RootsAnswer").innerHTML = "No real roots";
    }
    else if (mode=="temp"){
        return Func;
    }
}

function InitTrig(F,x,y,t)
{
  let tempLeft=Math.min(-1,Math.min(1.1 * t,0.91*t), Math.min(1.1 * x,0.91*x), Math.min(1.1 * y,0.91*y));
  let tempRight=Math.max(1,Math.max(1.1 * t,0.91*t), Math.max(1.1 * x,0.91*x), Math.max(1.1 * y,0.91*y));
  F.remapDomain(tempLeft,tempRight);
}

function DrawGraph(F, A, x, y, t) {
  canv = document.getElementById("graph");
  let ctx = canv.getContext("2d");
  let hgh = canv.height;
  let wth = canv.width;
  ctx.clearRect(0, 0, wth, hgh);
  let leftx;
  let rightx;
  let tempLeft=Math.min(-1,Math.min(1.1 * t,0.91*t), Math.min(1.1 * x,0.91*x), Math.min(1.1 * y,0.91*y));
  let tempRight=Math.max(1,Math.max(1.1 * t,0.91*t), Math.max(1.1 * x,0.91*x), Math.max(1.1 * y,0.91*y));
  let leftRoot;
  let rightRoot;
  if(F.rootsNum>20)
  {
      let i=0;
      while(F.roots[i+1]<tempLeft && typeof F.roots[i+1] !== undefined)
      {
          i++;
      }
      if(typeof F.roots[i+1]===undefined)
      {
          leftRoot=F.roots[i-1];
          rightRoot=F.roots[i];
      }
      else
      {
          if(F.roots[i+1]<tempRight){
              leftRoot=F.roots[i+1];
              i++;
          }
          else
              leftRoot=F.roots[i];
          while(F.roots[i+1]<tempRight && typeof F.roots[i+1]!==undefined)
              i++;
          if(typeof F.roots[i+1]===undefined)
              rightRoot=F.roots[i];
          else if(F.roots[i]>leftRoot)
              rightRoot=F.roots[i];
          else
              rightRoot=F.roots[i+1];
      }
  }
  else if(F.rootsNum>0)
  {
      leftRoot=F.roots[0];
      rightRoot=F.roots[F.rootsNum-1];
  }
  if(F.rootsNum>0){
    leftx = Math.min(Math.min(1.1 * leftRoot,0.91*leftRoot), tempLeft);
    rightx = Math.max(Math.max(1.1 * rightRoot,0.91*rightRoot), tempRight);
  }
  else {
    leftx=Math.min(-1, Math.min(1.1 * t,0.91*t), Math.min(1.1 * x,0.91*x), Math.min(1.1 * y,0.91*y));
    rightx = Math.max(1, Math.max(1.1 * t,0.91*t), Math.max(1.1 * x,0.91*x), Math.max(1.1 * y,0.91*y));
  }
  let fM=F.FindMaximum(leftx, rightx);
  let topy = Math.max(1, Math.max(1.1 * fM,0.91*fM));
  fM=F.FindMinimum(leftx, rightx);
  let bottomy = Math.min(-1, Math.min(1.1 * fM,0.91*fM));
  ctx.strokeStyle = "rgba(238,255,188,0.7)"; //#EEFFBC
  plotIntegral(F, x, y);
  ctx.beginPath();
  ctx.strokeStyle = "#000000";
  ctx.moveTo(0, scaleY(0));
  ctx.lineTo(wth, scaleY(0));
  ctx.stroke();
  ctx.closePath();
  ctx.beginPath();
  ctx.moveTo(scaleX(0), 0);
  ctx.lineTo(scaleX(0), hgh);
  ctx.stroke();
  ctx.closePath();
    
  ctx.font="12px Arial";
  let xMarkPlace=scaleY(0)-5;
  let yMarkPlace=scaleX(0)-5;
  ctx.moveTo(0.08*wth,xMarkPlace);
  ctx.lineTo(0.08*wth,xMarkPlace+10);
  ctx.stroke();
  ctx.fillText(antiScaleX(0.08*wth),0.08*wth-12,xMarkPlace-5);
  ctx.moveTo(0.92*wth,xMarkPlace);
  ctx.lineTo(0.92*wth,xMarkPlace+10);
  ctx.stroke();
  ctx.fillText(antiScaleX(0.92*wth),0.92*wth-12,xMarkPlace-5);
  ctx.moveTo(yMarkPlace,0.06*hgh);
  ctx.lineTo(yMarkPlace+10,0.06*hgh);
  ctx.stroke();
  ctx.fillText(antiScaleY(0.06*hgh),yMarkPlace+15,0.06*hgh+3);
  ctx.moveTo(yMarkPlace,0.94*hgh);
  ctx.lineTo(yMarkPlace+10,0.94*hgh);
  ctx.stroke();
  ctx.fillText(antiScaleY(0.94*hgh),yMarkPlace+15,0.94*hgh+3);
  
  ctx.strokeStyle = "#005169";
  plotFunction(F);
  ctx.strokeStyle = "rgb(113,35,76)";
  plotFunction(A);
    
  F.roots.forEach(element => {if(element>=leftRoot && element <= rightRoot) F.visibleRoots.push(element);});

  function plotFunction(G) {
    ctx.moveTo(scaleX(leftx), scaleY(G.Val(leftx)));
    ctx.beginPath();
    let dx = (rightx - leftx) / 10000;
    for (let i = leftx + dx; i <= rightx; i += dx) {
      if(G.isInDomain(i))
        ctx.lineTo(scaleX(i), scaleY(G.Val(i)));
      else
      {
        ctx.stroke();
        ctx.closePath();
        ctx.beginPath();
      }
    }
    ctx.stroke();
    ctx.closePath();
  }
  function plotIntegral(G, a, b) {
    let dx = (rightx - leftx) / 1000;
    let curFunc;
    ctx.beginPath();
    for (let i = a; i <= b; i += dx) {
      ctx.moveTo(scaleX(i), scaleY(0));
      curFunc = G.Val(i);
      ctx.lineTo(scaleX(i), scaleY(curFunc));
      ctx.stroke();
    }
    ctx.closePath();
  }
  function scaleX(n) {
    return (n - leftx) * wth / (rightx - leftx);
  }
  function scaleY(n) {
    return (topy - n) * hgh / (topy - bottomy);
  }
  function antiScaleX(n) {
    return n*(rightx-leftx)/wth+leftx;
  }
  function antiScaleY(n) {
    return topy-n*(topy-bottomy)/hgh;
  }
}

function drawMandelbrot(leftX, rightX, bottomY, topY, prec, itSelect, coloring) {
  let maxIt = 40;
  let precSelection;
  switch (prec) {
    case "1": precSelection = 300; maxIt = 30; break;
    case "2": precSelection = 500; maxIt = 40; break;
    case "3": precSelection = 650; maxIt = 60; break;
    case "4": precSelection = 800; maxIt = 100; break;
    case "5": precSelection = 1000; maxIt = 150; break;
  }
  switch (itSelect) {
    case "1": maxIt *= 0.5; break;
    case "2": break;
    case "3": maxIt *= 2; break;
    case "4": maxIt *= 5; break;
    case "5": maxIt *= 15; break;
  }
  canv = document.getElementById("mandelbrot");
  let ctx = canv.getContext("2d");
  let hgh = canv.height;
  let wth = canv.width;
  ctx.clearRect(0, 0, wth, hgh);
  let dx = (rightX - leftX) / precSelection;
  let dy = (bottomY - topY) / precSelection;
  let color = new Array(maxIt + 1);
  let c;
  let colorSystem;
  let colorFormat;
  let scaleddx = scaleX(dx);
  let scaleddy = scaleY(dy);
  let dtx = Math.max(1, scaleddx);
  let dty = Math.max(1, scaleddy);
    
  if(document.getElementById("mandelbrotDownloadButton")!=null)
      mandelbrotDownloadButton.remove();

  switch (coloring) {
    case "1":
      //Greyscale
      colorSystem="rgb(";
      colorFormat="";
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * i / maxIt;
        color[i].g = 255 * i / maxIt;
        color[i].b = 255 * i / maxIt;
      }
      break;
    case "2":
      colorSystem="rgb(";
      colorFormat="";
      //purplescale
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * i / maxIt;
        color[i].g = 255 * i / maxIt;
        color[i].b = 255 * Math.log(i) / Math.log(maxIt);
      }
      break;
    case "3":
      colorSystem="rgb(";
      colorFormat="";
      //red-background/shallow g,b
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * Math.log(i) / Math.log(maxIt);
        color[i].g = Math.max(0, -i * i + 2 * Math.sqrt(255) * i);
        color[i].b = Math.max(0, -2 * i * i + 2 * Math.sqrt(510) * i);
      }
      break;
    /*for(let i=0;i<=maxIt;i++)
        {
            color[i]=new Object();
            color[i].r=-1*Math.abs(Math.abs(x)-127.5)+127.5;
            color[i].g=Math.max(0,-i*i+2*Math.sqrt(255)*i);
            color[i].b=Math.max(0,-2*i*i+2*Math.sqrt(510)*i);
        }*/
    case "4":
      colorSystem="hsl(";
      colorFormat="%";
      //HSL purple-reddish
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 360 * Math.log(i) / Math.log(maxIt);
        color[i].g = 100 * i / maxIt;
        color[i].b = 75 * i / maxIt;
      }
      break;
    case "5":
      colorSystem="hsl(";
      colorFormat="%";
      //HSL normal
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 360 * i / maxIt;
        color[i].g = 100 * i / maxIt;
        color[i].b = 75 * i / maxIt;
      }
      break;
    case "6":
      colorSystem="rgb(";
      colorFormat="";
      //RGB sinusoid
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * 0.5*Math.sin(Math.PI*i/maxIt)+1;
        color[i].g = 255 * 0.5*Math.cos(Math.PI*i/maxIt)+1;
        color[i].b = 255 * Math.tan(Math.PI/4*i/maxIt);
      }
      break;
    case "7":
      colorSystem="rgb(";
      colorFormat="";
      //RGB sinusoid II
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * 0.5*Math.sin(10*Math.PI*i/maxIt)+1;
        color[i].g = 255 * 0.5*Math.cos(10*Math.PI*i/maxIt)+1;
        color[i].b = 255 * Math.tan(Math.PI/4*i/maxIt);
      }
      break;
    case "8":
      colorSystem="hsl(";
      colorFormat="%";
      //HSL sinusoid
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 360 * 0.5*Math.sin(10*Math.PI*i/maxIt)+1;
        color[i].g = 100 * 0.5*Math.cos(10*Math.PI*i/maxIt)+1;
        color[i].b = 75 * Math.tan(Math.PI/4*i/maxIt);
      }
      break;
  }

  for (let i = leftX; i <= rightX; i += dx)
    for (let j = topY; j >= bottomY; j += dy) {
      Mandelbrot(i, j, maxIt/*,i,j*/);
    }
    
    setupDownloadImage(canv,"mandelbrot");
    
  function Mandelbrot(x, y, i) {
    let a = x;
    let b = y;
    let temp;
    //let c;
    while (i > 0) {
      if (a * a + b * b > 4) {
        c = colorSystem + color[i].r + ',' + color[i].g + colorFormat+',' + color[i].b + colorFormat+')';
        ctx.fillStyle = c;
        ctx.fillRect(scaleX(x), scaleY(y), dtx, dty);
        break;
      }
      else {
        temp = a;
        a = a * a - b * b + x;
        b = 2 * temp * b + y;
        i--;
      }
    }
    if (i <= 0) {
      ctx.fillStyle = "rgb(0,0,0)";
      ctx.fillRect(scaleX(x), scaleY(y), dtx, dty);
    }
  }
  function scaleX(n) {
    return (n - leftX) * wth / (rightX - leftX);
  }
  function scaleY(n) {
    return (topY - n) * hgh / (topY - bottomY);
  }
}

function lockProportions(c1,c2,c3,mode)
{
  if(mode==2)
    return (c3-c2)+c1;
  else if(mode==1)
    return c1-(c3-c2);
  else if(mode==3)
    return c3-(c2-c1);
  else //mode==4
    return c3+(c2-c1);
}

function setupDownloadImage(canv,idFunc)
{
  let ctx = canv.getContext("2d");
  let hgh = canv.height;
  let wth = canv.width;
  let img = canv.toDataURL();
  document.createElement("br");
  var button = document.createElement("input"); //Create <a>
  button.type="button";
  button.id=idFunc+"DownloadButton";
  //canv.insertAdjacentElement("afterend", button);
  button.value="Download as PNG";
  button.addEventListener('click', downloadImage.bind(this, img,idFunc));
  //button.onclick="downloadImage(img)";
  canv.insertAdjacentElement("afterend", button);
  //button.onclick = "downloadImage(img,\"mandelbrot\")"; //Image Base64 Goes here
  //a.download = "Image.png"; //File name Here
  //a.click();;
}

function downloadImage(img,filename)
{
    //let i=img;
    var a = document.createElement("a"); //Create <a>
    a.href = img; //Image Base64 Goes here
    a.download = filename+".png"; //File name Here
    //a.download="mandelbrot.png";
    a.click();
}

//julia

function drawJulia(leftX, rightX, bottomY, topY, pa,pb, prec, itSelect, coloring) {
  let maxIt = 40;
  let precSelection;
  switch (prec) {
    case "1": precSelection = 300; maxIt = 30; break;
    case "2": precSelection = 500; maxIt = 40; break;
    case "3": precSelection = 650; maxIt = 60; break;
    case "4": precSelection = 800; maxIt = 100; break;
    case "5": precSelection = 1000; maxIt = 150; break;
  }
  switch (itSelect) {
    case "1": maxIt *= 0.5; break;
    case "2": break;
    case "3": maxIt *= 2; break;
    case "4": maxIt *= 5; break;
    case "5": maxIt *= 15; break;
  }
  canv = document.getElementById("julia");
  let ctx = canv.getContext("2d");
  let hgh = canv.height;
  let wth = canv.width;
  ctx.clearRect(0, 0, wth, hgh);
  let dx = (rightX - leftX) / precSelection;
  let dy = (bottomY - topY) / precSelection;
  let color = new Array(maxIt + 1);
  let c;
  let colorSystem;
  let colorFormat;
  let scaleddx = scaleX(dx);
  let scaleddy = scaleY(dy);
  let dtx = Math.max(1, scaleddx);
  let dty = Math.max(1, scaleddy);
    
  if(document.getElementById("juliaDownloadButton")!=null)
      juliaDownloadButton.remove();

  switch (coloring) {
    case "1":
      //Greyscale
      colorSystem="rgb(";
      colorFormat="";
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * i / maxIt;
        color[i].g = 255 * i / maxIt;
        color[i].b = 255 * i / maxIt;
      }
      break;
    case "2":
      colorSystem="rgb(";
      colorFormat="";
      //purplescale
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * i / maxIt;
        color[i].g = 255 * i / maxIt;
        color[i].b = 255 * Math.log(i) / Math.log(maxIt);
      }
      break;
    case "3":
      colorSystem="rgb(";
      colorFormat="";
      //red-background/shallow g,b
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * Math.log(i) / Math.log(maxIt);
        color[i].g = Math.max(0, -i * i + 2 * Math.sqrt(255) * i);
        color[i].b = Math.max(0, -2 * i * i + 2 * Math.sqrt(510) * i);
      }
      break;
    /*for(let i=0;i<=maxIt;i++)
        {
            color[i]=new Object();
            color[i].r=-1*Math.abs(Math.abs(x)-127.5)+127.5;
            color[i].g=Math.max(0,-i*i+2*Math.sqrt(255)*i);
            color[i].b=Math.max(0,-2*i*i+2*Math.sqrt(510)*i);
        }*/
    case "4":
      colorSystem="hsl(";
      colorFormat="%";
      //HSL purple-reddish
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 360 * Math.log(i) / Math.log(maxIt);
        color[i].g = 100 * i / maxIt;
        color[i].b = 75 * i / maxIt;
      }
      break;
    case "5":
      colorSystem="hsl(";
      colorFormat="%";
      //HSL normal
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 360 * i / maxIt;
        color[i].g = 100 * i / maxIt;
        color[i].b = 75 * i / maxIt;
      }
      break;
    case "6":
      colorSystem="rgb(";
      colorFormat="";
      //RGB sinusoid
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * 0.5*Math.sin(Math.PI*i/maxIt)+1;
        color[i].g = 255 * 0.5*Math.cos(Math.PI*i/maxIt)+1;
        color[i].b = 255 * Math.tan(Math.PI/4*i/maxIt);
      }
      break;
    case "7":
      colorSystem="rgb(";
      colorFormat="";
      //RGB sinusoid II
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 255 * 0.5*Math.sin(10*Math.PI*i/maxIt)+1;
        color[i].g = 255 * 0.5*Math.cos(10*Math.PI*i/maxIt)+1;
        color[i].b = 255 * Math.tan(Math.PI/4*i/maxIt);
      }
      break;
    case "8":
      colorSystem="hsl(";
      colorFormat="%";
      //HSL sinusoid
      for (let i = 0; i <= maxIt; i++) {
        color[i] = new Object();
        color[i].r = 360 * 0.5*Math.sin(10*Math.PI*i/maxIt)+1;
        color[i].g = 100 * 0.5*Math.cos(10*Math.PI*i/maxIt)+1;
        color[i].b = 75 * Math.tan(Math.PI/4*i/maxIt);
      }
      break;
  }

  for (let i = leftX; i <= rightX; i += dx)
    for (let j = topY; j >= bottomY; j += dy) {
      Julia(i, j, maxIt,pa,pb);
    }
    
    setupDownloadImage(canv,"julia");
    
  function Julia(x, y, i,pa,pb) {
    let a = x;
    let b = y;
    let temp;
    //let c;
    while (i > 0) {
      if (a * a + b * b > 4) {
        c = colorSystem + color[i].r + ',' + color[i].g + colorFormat+',' + color[i].b + colorFormat+')';
        ctx.fillStyle = c;
        ctx.fillRect(scaleX(x), scaleY(y), dtx, dty);
        break;
      }
      else {
        temp = a;
        a = a * a - b * b + pa;
        b = 2 * temp * b + pb;
        i--;
      }
    }
    if (i <= 0) {
      ctx.fillStyle = "rgb(0,0,0)";
      ctx.fillRect(scaleX(x), scaleY(y), dtx, dty);
    }
  }
  function scaleX(n) {
    return (n - leftX) * wth / (rightX - leftX);
  }
  function scaleY(n) {
    return (topY - n) * hgh / (topY - bottomY);
  }
}

function selectExample(setting,scope){
    let exCoordinates = new Array(6);
    exCoordinates=setting.split(',').map(Number);
    
    document.getElementById(scope+"MinX").value=exCoordinates[0];
    document.getElementById(scope+"MaxX").value=exCoordinates[1];
    document.getElementById(scope+"MinY").value=exCoordinates[2];
    document.getElementById(scope+"MaxY").value=exCoordinates[3]; 
    
    if(document.getElementById(scope+"ParameterA")!=null)
        document.getElementById(scope+"ParameterA").value=exCoordinates[4];
    
    if(document.getElementById(scope+"ParameterB")!=null)
        document.getElementById(scope+"ParameterB").value=exCoordinates[5];
}