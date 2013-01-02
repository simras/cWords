
class Array
  
  # threaded each
  def threach(n = 1, &b)
    return [] if n == 0 or size == 0
    result = Array.new(size)
    return self.send(:each,&b) if n == 1 # trying return here
    
    n = [n,size].min
    
    part_size, part_remainder = size/n, size % n
    threads = []

    pstart = 0
    n.times do |pi|
      pend = pstart + part_size - 1
      pend += 1 if pi<part_remainder
      threads << Thread.new(pstart,pend) do |a,b|
        for j in a..b
          yield(slice(j))
        end
      end
      pstart = pend+1
    end
    
    threads.each { |t| t.join }
    self
  end
  
  # unit tests for threach
  #Array.send(:include, Threach)
  #a=(0..4).to_a
  #res = a.map{|x| x*10}
  #a.size.times do |pn|
  #  b=Array.new(a.size)
  #  a.threach(pn+1) {|x| b[x]=x*10}
  #  puts b == beach
  #end

  
  def shuffle()
    arr = self.dup
    arr.size.downto 2 do |j|
      r = Kernel::rand(j)
      arr[j-1], arr[r] = arr[r], arr[j-1]
    end
    arr
  end
  
  #reference: http://blade.nagaokaut.ac.jp/~sinara/ruby/math/combinatorics/array-rep_perm.rb
  def rep_perm(n)
    if n < 0
    elsif n == 0
      yield([])
    else
      rep_perm(n - 1) do |x|
        each do |y|
            yield(x + [y])
          end
      end
    end
  end

  def to_statarray
    StatArray.new(self)
  end

  #return the (sorted) ranks of the elements
  def ranks()
    h = Hash.new {|h,k| h[k]=[]}
    self.sort.each_with_index{|x,idx| h[x] << idx}
    self.map{|x| h[x].first + (h[x].size-1)/2.0}
  end
  
end

class StatArray < Array

  alias :count size

  def sum
    inject(0) { |sum, x| sum + x }
  end

  def mean
    return 0.0 if self.size == 0
    sum.to_f / self.size
  end
  alias :arithmetic_mean :mean

  def median
    return 0 if self.size == 0
    tmp = sort
    mid = tmp.size / 2
    if (tmp.size % 2) == 0
      (tmp[mid-1] + tmp[mid]).to_f / 2
    else
      tmp[mid]
    end
  end

  # The sum of the squared deviations from the mean.
  def summed_sqdevs
    return 0 if count < 2
    m = mean
    StatArray.new(map { |x| (x - m) ** 2 }).sum
  end

  # Variance of the sample.
  def variance
    # Variance of 0 or 1 elements is 0.0
    return 0.0 if count < 2
    summed_sqdevs / (count - 1)
  end

  # Variance of a population.
  def pvariance
    # Variance of 0 or 1 elements is 0.0
    return 0.0 if count < 2
    summed_sqdevs / count
  end

  # Standard deviation of a sample.
  def stddev
    Math::sqrt(variance)
  end

  # Standard deviation of a population.
  def pstddev
    Math::sqrt(pvariance)
  end

  # Calculates the standard error of this sample.
  def stderr
    return 0.0 if count < 2
    stddev/Math::sqrt(size)
  end

  def gmean
    return 0.0 if self.size == 0
    return nil if self.any?{|x| x < 0 } # not negative
    return 0.0 if self.any?{|x| x == 0 } #includes 0
    return 10**(self.map{|x| log10(x)}.sum/self.size)
  end

  #rank cdf - array is not sorted before cdf calculation
  def rcdf(normalised = 1.0)
    s = sum.to_f
    inject([0.0]) { |c,d| c << c[-1] + normalised.to_f*d.to_f/s }
  end

end


class String
  
  def shuffle
    self.split("").shuffle.join
  end

  def to_f2
    return 1/0.0 if self == "inf"
    return -1/0.0 if self == "-inf"
    return 0/0.0 if self == "NaN"
    self.to_f
  end
  
end


class Float
  
  alias_method :orig_to_s, :to_s

  def to_s(arg = nil)
    if arg.nil?
      #orig_to_s
      sprintf("%f", self)
    else
      sprintf("%.#{arg}f", self)
    end
  end

  def to_e(arg = nil)
    if arg.nil?
      #orig_to_s
      sprintf("%.#{2}e", self)
    else
      sprintf("%.#{arg}e", self)
    end
  end

  def sign
    self < 0 ? -1 : 1
  end
  
end

def log2(number)
  Math.log(number)/Math.log(2)
end

def log10(number)
  Math.log(number)/Math.log(10)
end

def N2num(n)
  if n == 'A' then return 0 end
  if n == 'C' then return 2 end
  if n == 'G' then return 1 end
  if n == 'T' then return 3
  else 
    print "Error N2num wordRS ","\n"
    return -1
  end
end


def num2N(n)
  case n
  when 0 then return 'A'
  when 2 then return 'C'
  when 1 then return 'G'
  when 3 then return 'T'
  else 
    print "Error num2N wordRS","\n"
    return 'X' 
  end
end

def calc_idx(s)
  ss = s.split("")
  k = s.size
  sum = 0
  if k == 1 then
    return N2num(ss[0])
  end
  for i in 0..(k-1) do
    sum = sum + (4**(k - i - 1)) * N2num(ss[i])
  end
  return sum
end

def calc_patt2(idx ,wss)
  j = 0
  summ = 4**wss[j]
  while idx >= summ do
    j = j + 1
    summ = summ + 4**wss[j]
  end
  ws = wss[j]
  old_summ = summ - 4**ws
  return calc_patt(idx - old_summ, ws)
end

def calc_patt(n, k)
    st = ""
  sum = 0;
  for i in 0..(k-1) do
      temp = (n / 4**(k - i - 1))
    n = n - temp * 4**(k - i - 1)
      st = st + num2N(temp)
  end
  return st
end

def calcMean(sig,t)
  return (Math.sqrt(Math.PI)*sig*Math.log(2))/(Math.sqrt(2.0/t))
  end

def calcVar(sig,t,m)
  return t * sig**2 * ((Math.PI**2)/12) - m**2
end
  
def calc_pval(x,sig,t)
  sum=1
  if(x==0) then
    return 1;
  end
  for h in 1..10000 do
    p = (-1**h)*Math.exp((-2*(x**2)*(h**2))/((sig**2)*t))
    sum = sum + 2*p
    if (p.abs < 10**-10) then
      if (1-sum==0) then
        return 10**-10
      else
        return 1-sum
      end
    end
  end
  print "Error: Variance did not converge" , "\n"
  if 1-sum==0 then
    return (10**-10)
  else
    return (1-sum)
  end
end

def random_string(length=5)
  chars = 'abcdefghjkmnpqrstuvwxyzABCDEFGHJKLMNPQRSTUVWXYZ23456789'
  password = ''
  length.times { password << chars[rand(chars.size)] }
  return password
end


def calc_idx_ext(patt,ws)
  done = false
  sum = 0
  i = 0
  k = patt.size
  while !done
    if k == ws[i]
      done = true
    else
      sum = sum + 4**ws[i]
      i = i + 1
    end
  end
  return sum + calc_idx(patt)
end

def calcSig(x)
  # This function calculates the squared sigma
  # which is a parameter of the Rayleigh Distribution
  n = x.length
  sum = 0
  x.each do |val|
    sum = sum + val**2
  end
  return sum/(2.0*n)
end

def rayleighMean(sig)
  # This function calculates the Mean of the Rayleigh Distribution
  return Math.sqrt((sig*Math::PI)/2.0)
end
    
def rayleighVar(sig)
  # This function calculates the Variance of the Rayleigh Distribution
  return ((4.0 - Math::PI) * sig)/2.0
end

def rayleighPval(x,sig)
  # This function calculates the P-vals of the Rayleigh Distribution
  return Math.exp(-(x**2/(sig)))
end

def cmp(a,b,i)
  ## Compares two 2d array elements by 1st element
  return 1  if a[i] < b[i]
  return 0  if a[i] == b[i]
  return -1 if a[i] > b[i]
end

def minimize(obsVals,funArray,n,m,test)
  # n: number of rows
  # m: number of columns
  # obsVals: the observed values we will fit to the exponential running sum
  # funArray: precalculated theoretical eksponential running sum

 # print "Length1 ", funArray.size(), " Length2 ",funArray[1499].size() ,"\n"
  minIdx = 0
  minimumSqError = 10000000
  minArray = Array.new(20)
  for i in (0..(n-1))
    sqerror = 0
    for j in (0..(m-1))
      error = obsVals[j] - funArray[i][j]
      sqerror = sqerror + error * error
    end
 #   print "min Indx ",sqerror,"\n"
    if minimumSqError > sqerror then
      minimumSqError = sqerror
      minIdx = i
      if "test" == test then
        minArray = funArray[i]
      end
#      print "min Indx ", minIdx, " minSQErr ", minimumSqError, "\n"
    end
  end
  if "test" == test then
    print "Test ",minIdx," \n", obsVals.join(" "),"\n"
    print minArray.join(" ") ,"\n"
  end
  return [minIdx,minimumSqError]
end

def normalize(valArray,maxrs)
  for i in (0..(valArray.length-1))
    valArray[i] = (valArray[i])/(maxrs)
  end
  return valArray
end

def chooseValues(valArray,maxrs)
  # Select values in value array
  # with even spacing
  # 
  #normalize(valArray,pmean,z_sc)
  reducArr = Array.new(19)
  n = valArray.length()
  step = (n/20).floor
  for i in (1..19)
   # print maxrs,"\n"
    if (i*step)-1 > n then
      # normalize: div by maxrs
      reducArr[i-1] =  valArray[n-1]/maxrs
    else
      reducArr[i-1] = valArray[(i * step)-1]/maxrs
    end
  end
#  print reducArr.length()
  return reducArr
end

def fitFun(valArray,maxrs,test)
  # 
  # 
  # 
  
  obsVals = chooseValues(valArray,maxrs)
  alphaStep = 0.01
  n=150
  m=20
  # open fun file
  funFile = File.open("resources/fitFun.txt")
  funArr = Array.new(2*n){Array.new(m-1)}
  i = 0
  lineList = []
  while (line = funFile.gets)
    lineList = line.split()
    for k in (0..18)
      funArr[i][k] = lineList[k].to_f
    end
#    print funArr[i],"\n"
    i = i + 1
  end
  n = i
  m = lineList.length()
  # close fun file
  funFile.close()
 # print "Number of values ",m,"\n"
  minimizeList = minimize(obsVals,funArr,n,m,test)
  alphaNum = minimizeList[0]
  sumSqError = minimizeList[1]
  if alphaNum < 150 then
    alpha = 1 + alphaStep + alphaNum * alphaStep
  else
    alpha = -(1 + alphaStep + (alphaNum-150) * alphaStep)
  end
 # print "Alpha ",alpha," Sum of Squared errors ",sumSqError,"\n"
  
  return [alpha,sumSqError]
end
