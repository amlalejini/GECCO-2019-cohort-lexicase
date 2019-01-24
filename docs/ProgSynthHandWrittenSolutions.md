# Hand-coded solutions for programming synthesis benchmark problems

<!-- TOC -->

- [Problem - Number IO](#problem---number-io)
- [Problem - Small or Large](#problem---small-or-large)
- [Problem - ForLoopIndex](#problem---forloopindex)
- [Problem - CompareStringLengths](#problem---comparestringlengths)
- [Problem - CollatzNumbers](#problem---collatznumbers)
- [Problem - StringLengthsBackwards](#problem---stringlengthsbackwards)
- [Problem - LastIndexOfZero](#problem---lastindexofzero)
- [Problem - CountOdds](#problem---countodds)
- [Problem - Mirror Image](#problem---mirror-image)
- [Problem - Vectors Summed](#problem---vectors-summed)
- [Problem - Sum of Squares](#problem---sum-of-squares)
- [Problem - Vector Average](#problem---vector-average)
- [Problem - Median](#problem---median)
- [Problem - Smallest](#problem---smallest)
- [Problem - Grade](#problem---grade)

<!-- /TOC -->

## Problem - Number IO

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  sol.PushInst("LoadInt", {matrix[0], matrix[0], matrix[0]});
  sol.PushInst("LoadDouble", {matrix[1], matrix[1], matrix[1]});
  sol.PushInst("Add", {matrix[0], matrix[1], matrix[2]});
  sol.PushInst("SubmitNum", {matrix[2], matrix[2], matrix[2]});
  sol.PushInst("Return", {matrix[0], matrix[0], matrix[0]});
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - Small or Large

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  // Create thresholds to compare n to.
  sol.PushInst("Set-10", {matrix[0], matrix[4], matrix[4]});
  sol.PushInst("Set-2",  {matrix[1], matrix[4], matrix[4]});
  sol.PushInst("Mult",   {matrix[0], matrix[0], matrix[2]});
  sol.PushInst("Mult",   {matrix[0], matrix[2], matrix[0]});
  sol.PushInst("Mult",   {matrix[0], matrix[1], matrix[1]});
  // Load input
  sol.PushInst("LoadInt",     {matrix[2], matrix[4], matrix[4]});
  // Check if n < 1000
  sol.PushInst("TestNumLess", {matrix[2], matrix[0], matrix[3]});
  sol.PushInst("If",          {matrix[3], matrix[4], matrix[4]});
  sol.PushInst("SubmitSmall", {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Return",      {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Close",       {matrix[4], matrix[4], matrix[4]});
  // Check if n < 2000
  sol.PushInst("TestNumLess", {matrix[2], matrix[1], matrix[3]});
  sol.PushInst("If",          {matrix[3], matrix[4], matrix[4]});
  sol.PushInst("SubmitNone",  {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Return",      {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Close",       {matrix[4], matrix[4], matrix[4]});
  // n must be >= 2000
  sol.PushInst("SubmitLarge", {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Return",      {matrix[4], matrix[4], matrix[4]});
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - ForLoopIndex

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("CopyMem",     {matrix[0], matrix[4], matrix[7]});
  sol.PushInst("Inc",         {matrix[5], matrix[7], matrix[7]});
  sol.PushInst("While",       {matrix[5], matrix[7], matrix[7]});
  sol.PushInst("SubmitNum",   {matrix[4], matrix[7], matrix[7]});
  sol.PushInst("Add",         {matrix[4], matrix[2], matrix[4]});
  sol.PushInst("TestNumLess", {matrix[4], matrix[1], matrix[5]});
  sol.PushInst("Close",       {matrix[7], matrix[7], matrix[7]});
  sol.PushInst("Return",      {matrix[7], matrix[7], matrix[7]});

  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - CompareStringLengths

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("StrLength",   {matrix[0], matrix[0], matrix[4]});
  sol.PushInst("StrLength",   {matrix[1], matrix[1], matrix[4]});
  sol.PushInst("StrLength",   {matrix[2], matrix[2], matrix[4]});
  sol.PushInst("TestNumLess", {matrix[0], matrix[1], matrix[3]});
  sol.PushInst("If",          {matrix[3], matrix[4], matrix[4]});
  sol.PushInst("TestNumLess", {matrix[1], matrix[2], matrix[3]});
  sol.PushInst("If",          {matrix[3], matrix[4], matrix[4]});
  sol.PushInst("SubmitVal",   {matrix[3], matrix[4], matrix[4]});
  sol.PushInst("Close",       {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("Close",       {matrix[4], matrix[4], matrix[4]});
  sol.PushInst("SubmitVal",   {matrix[3], matrix[4], matrix[4]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - CollatzNumbers

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("Inc",           {matrix[1], matrix[8], matrix[8]});
  sol.PushInst("Set-2",         {matrix[2], matrix[8], matrix[8]});
  sol.PushInst("Set-3",         {matrix[3], matrix[8], matrix[8]});
  sol.PushInst("Inc",           {matrix[4], matrix[8], matrix[8]});
  sol.PushInst("TestNumNEqu",   {matrix[1], matrix[0], matrix[5]});
  sol.PushInst("While",         {matrix[5], matrix[8], matrix[8]});
  sol.PushInst(  "Inc",         {matrix[4], matrix[8], matrix[8]});
  sol.PushInst(  "Mod",         {matrix[0], matrix[2], matrix[6]});
  sol.PushInst(  "IfNot",       {matrix[6], matrix[8], matrix[8]});
  sol.PushInst(    "Div",       {matrix[0], matrix[2], matrix[0]});
  sol.PushInst(  "Close",       {matrix[8], matrix[8], matrix[8]});
  sol.PushInst(  "If",          {matrix[6], matrix[8], matrix[8]});
  sol.PushInst(    "Mult",      {matrix[0], matrix[3], matrix[0]});
  sol.PushInst(    "Inc",       {matrix[0], matrix[8], matrix[8]});
  sol.PushInst(  "Close",       {matrix[8], matrix[8], matrix[8]});
  sol.PushInst(  "TestNumNEqu", {matrix[1], matrix[0], matrix[5]});
  sol.PushInst("Close",         {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("SubmitNum",     {matrix[4], matrix[8], matrix[8]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - StringLengthsBackwards

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("LoadStrVec",  {matrix[0], matrix[4], matrix[4]});
  sol.PushInst("VecReverse",  {matrix[0], matrix[4], matrix[4]});
  sol.PushInst("Foreach",     {matrix[1], matrix[0], matrix[4]});
  sol.PushInst(  "StrLength", {matrix[1], matrix[2], matrix[4]});
  sol.PushInst(  "SubmitVal", {matrix[2], matrix[4], matrix[4]});
  sol.PushInst("Close",       {matrix[4], matrix[4], matrix[4]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - LastIndexOfZero

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("LoadVec",     {matrix[0], matrix[4], matrix[4]});
  sol.PushInst("VecReverse",  {matrix[0], matrix[4], matrix[4]});
  sol.PushInst("VecIndexOf",  {matrix[0], matrix[1], matrix[2]});
  sol.PushInst("VecLen",      {matrix[0], matrix[1], matrix[4]});
  sol.PushInst("Sub",         {matrix[1], matrix[2], matrix[3]});
  sol.PushInst("Dec",         {matrix[3], matrix[4], matrix[4]});
  sol.PushInst("SubmitNum",   {matrix[3], matrix[4], matrix[4]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - CountOdds

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("LoadVec",   {matrix[0], matrix[7], matrix[7]});
  sol.PushInst("Set-2",     {matrix[2], matrix[7], matrix[7]});
  sol.PushInst("Foreach",   {matrix[1], matrix[0], matrix[7]});
  sol.PushInst(  "Mod",     {matrix[1], matrix[2], matrix[3]});
  sol.PushInst(  "If",      {matrix[3], matrix[7], matrix[7]});
  sol.PushInst(    "Inc",   {matrix[4], matrix[7], matrix[7]});
  sol.PushInst(  "Close",   {matrix[7], matrix[7], matrix[7]});
  sol.PushInst("Close",     {matrix[7], matrix[7], matrix[7]});
  sol.PushInst("SubmitNum", {matrix[4], matrix[7], matrix[7]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - Mirror Image

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("LoadVec1",   {matrix[0], matrix[7], matrix[7]});
  sol.PushInst("LoadVec2",   {matrix[1], matrix[7], matrix[7]});
  sol.PushInst("VecReverse", {matrix[1], matrix[7], matrix[7]});
  sol.PushInst("TestMemEqu", {matrix[0], matrix[1], matrix[2]});
  sol.PushInst("SubmitVal",  {matrix[2], matrix[7], matrix[7]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - Vectors Summed

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("LoadVec1",      {matrix[0], matrix[8], matrix[8]});
  sol.PushInst("LoadVec2",      {matrix[1], matrix[8], matrix[8]});
  sol.PushInst("VecLen",        {matrix[0], matrix[2], matrix[8]});
  sol.PushInst("Set-0",         {matrix[3], matrix[8], matrix[8]});
  sol.PushInst("TestNumLess",   {matrix[3], matrix[2], matrix[4]});
  sol.PushInst("While",         {matrix[4], matrix[8], matrix[8]});
  sol.PushInst(  "VecGet",      {matrix[0], matrix[3], matrix[5]});
  sol.PushInst(  "VecGet",      {matrix[1], matrix[3], matrix[6]});
  sol.PushInst(  "Add",         {matrix[5], matrix[6], matrix[7]});
  sol.PushInst(  "VecSet",      {matrix[0], matrix[3], matrix[7]});
  sol.PushInst(  "Inc",         {matrix[3], matrix[8], matrix[8]});
  sol.PushInst(  "TestNumLess", {matrix[3], matrix[2], matrix[4]});
  sol.PushInst("Close",         {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("SubmitVec",     {matrix[0], matrix[8], matrix[8]});

  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - Sum of Squares

Closed form solution:

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("LoadNum",    {matrix[0], matrix[8], matrix[8]});
  sol.PushInst("Set-2",      {matrix[1], matrix[8], matrix[8]});
  sol.PushInst("Set-6",      {matrix[2], matrix[8], matrix[8]});
  sol.PushInst("CopyMem",    {matrix[0], matrix[3], matrix[8]});
  sol.PushInst("Inc",        {matrix[3], matrix[8], matrix[8]});
  sol.PushInst("Mult",       {matrix[0], matrix[3], matrix[4]});
  sol.PushInst("Mult",       {matrix[0], matrix[1], matrix[5]});
  sol.PushInst("Inc",        {matrix[5], matrix[8], matrix[8]});
  sol.PushInst("Mult",       {matrix[4], matrix[5], matrix[6]});
  sol.PushInst("Div",        {matrix[6], matrix[2], matrix[7]});
  sol.PushInst("SubmitNum",  {matrix[7], matrix[8], matrix[8]});
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - Vector Average

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("LoadVec",   {matrix[0], matrix[5], matrix[5]});
  sol.PushInst("Foreach",   {matrix[1], matrix[0], matrix[5]});
  sol.PushInst(  "Add",     {matrix[1], matrix[2], matrix[2]});
  sol.PushInst("Close",     {matrix[5], matrix[5], matrix[5]});
  sol.PushInst("VecLen",    {matrix[0], matrix[3], matrix[5]});
  sol.PushInst("Div",       {matrix[2], matrix[3], matrix[4]});
  sol.PushInst("SubmitNum", {matrix[4], matrix[5], matrix[5]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - Median

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("MakeVector",    {matrix[0], matrix[2], matrix[4]});
  sol.PushInst("LoadNum1",      {matrix[5], matrix[0], matrix[0]});
  sol.PushInst("LoadNum1",      {matrix[6], matrix[0], matrix[0]});
  sol.PushInst("Foreach",       {matrix[0], matrix[4], matrix[0]});
  sol.PushInst(  "TestNumLess", {matrix[0], matrix[5], matrix[1]});
  sol.PushInst(  "If",          {matrix[1], matrix[0], matrix[0]});
  sol.PushInst(    "CopyMem",   {matrix[0], matrix[5], matrix[0]});
  sol.PushInst(  "Close",       {matrix[0], matrix[0], matrix[0]});
  sol.PushInst(  "TestNumLess", {matrix[6], matrix[0], matrix[1]});
  sol.PushInst(  "If",          {matrix[1], matrix[0], matrix[0]});
  sol.PushInst(    "CopyMem",   {matrix[0], matrix[6], matrix[0]});
  sol.PushInst(  "Close",       {matrix[0], matrix[0], matrix[0]});
  sol.PushInst(  "Add",         {matrix[0], matrix[7], matrix[7]});
  sol.PushInst("Close",         {matrix[0], matrix[0], matrix[0]});
  sol.PushInst("Sub",           {matrix[7], matrix[5], matrix[7]});
  sol.PushInst("Sub",           {matrix[7], matrix[6], matrix[7]});
  sol.PushInst("SubmitNum",     {matrix[7], matrix[0], matrix[0]});

  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - Smallest

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);

  sol.PushInst("LoadNum1",    {matrix[0], matrix[7], matrix[7]});
  sol.PushInst("LoadNum2",    {matrix[1], matrix[7], matrix[7]});
  sol.PushInst("LoadNum3",    {matrix[2], matrix[7], matrix[7]});
  sol.PushInst("LoadNum4",    {matrix[3], matrix[7], matrix[7]});
  sol.PushInst("MakeVector",  {matrix[0], matrix[3], matrix[4]});
  sol.PushInst("Foreach",     {matrix[5], matrix[4], matrix[7]});
  sol.PushInst("TestNumLess", {matrix[5], matrix[0], matrix[6]});
  sol.PushInst("If",          {matrix[6], matrix[7], matrix[7]});
  sol.PushInst("CopyMem",     {matrix[5], matrix[0], matrix[7]});
  sol.PushInst("Close",       {matrix[7], matrix[7], matrix[7]});
  sol.PushInst("Close",       {matrix[7], matrix[7], matrix[7]});
  sol.PushInst("SubmitNum",   {matrix[0], matrix[7], matrix[7]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```

## Problem - Grade

```{C++}
  emp::vector<emp::BitSet<TAG_WIDTH>> matrix = GenHadamardMatrix<TAG_WIDTH>();
  hardware_t::Program sol(inst_lib);
  sol.PushInst("LoadThreshA",        {matrix[0], matrix[8], matrix[8]});
  sol.PushInst("LoadThreshB",        {matrix[1], matrix[8], matrix[8]});
  sol.PushInst("LoadThreshC",        {matrix[2], matrix[8], matrix[8]});
  sol.PushInst("LoadThreshD",        {matrix[3], matrix[8], matrix[8]});
  sol.PushInst("LoadGrade",          {matrix[4], matrix[8], matrix[8]});
  sol.PushInst("TestNumGreaterTEqu", {matrix[4], matrix[0], matrix[5]});
  sol.PushInst("If",                 {matrix[5], matrix[8], matrix[8]});
    sol.PushInst("SubmitA",          {matrix[8], matrix[8], matrix[8]});
    sol.PushInst("Return",           {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("Close",              {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("TestNumGreaterTEqu", {matrix[4], matrix[1], matrix[5]});
  sol.PushInst("If",                 {matrix[5], matrix[8], matrix[8]});
    sol.PushInst("SubmitB",          {matrix[8], matrix[8], matrix[8]});
    sol.PushInst("Return",           {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("Close",              {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("TestNumGreaterTEqu", {matrix[4], matrix[2], matrix[5]});
  sol.PushInst("If",                 {matrix[5], matrix[8], matrix[8]});
    sol.PushInst("SubmitC",          {matrix[8], matrix[8], matrix[8]});
    sol.PushInst("Return",           {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("Close",              {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("TestNumGreaterTEqu", {matrix[4], matrix[3], matrix[5]});
  sol.PushInst("If",                 {matrix[5], matrix[8], matrix[8]});
    sol.PushInst("SubmitD",          {matrix[8], matrix[8], matrix[8]});
    sol.PushInst("Return",           {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("Close",              {matrix[8], matrix[8], matrix[8]});
  sol.PushInst("SubmitF",            {matrix[8], matrix[8], matrix[8]});
  
  prog_world->Inject(sol, PROG_POP_SIZE);
```