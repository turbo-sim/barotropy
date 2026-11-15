/* Run with node --test benchmarks/test_report.cjs; no npm install needed. */
const {test} = require('node:test');
const assert = require('node:assert/strict');
const report = require('./report.js');
const sample = (a,b,c,d,e) => ({construct:a,solve:b,fit:c,evaluate_501:d,export:e});
const measurement = {case:{fluid:'CO2'},parameters:{b:2,a:1},samples:[sample(5,5,5,5,5),sample(1,100,1,1,1),sample(1,1,100,1,1),sample(100,1,1,1,1)]};
test('total median is the median of complete repetitions, not sum of stage medians',()=>{
  assert.equal(report.stats(measurement,'total').median,104);
  assert.equal(report.stats(measurement,'total').first,25);
  assert.equal(report.stats(measurement,'solve').median,1);
  assert.equal(report.stats(measurement,'total').repeats,3);
});
test('missing and modified cases cannot produce baseline ratios',()=>{
  const a={cases:{case:measurement}}, b={cases:{case:{...measurement,parameters:{a:1,b:2}}}};
  assert.equal(report.compatible(a,b,'case'),true);
  b.cases.case.parameters.a=9;
  assert.equal(report.compatible(a,b,'case'),false);
  assert.equal(report.compatible(a,{cases:{}},'case'),false);
});
test('negative and nonfinite timings are rejected',()=>{
  const data={schema_version:1,cases:{case:structuredClone(measurement)}};
  assert.equal(report.validate(data),data);
  data.cases.case.samples[1].solve=NaN;
  assert.throws(()=>report.validate(data),/invalid solve/);
  data.cases.case.samples[1].solve=-1;
  assert.throws(()=>report.validate(data),/invalid solve/);
});
test('fit summaries ignore undefined relative errors at zero reference',()=>{
  assert.deepEqual(report.fitErrorSummary({absolute_error:[-.1,.2,null],relative_error_percent:[null,-5,10]}),
    {maximumAbsolute:.2,maximumRelative:10});
  assert.deepEqual(report.fitErrorSummary({absolute_error:[0,0],relative_error_percent:[null,null]}),
    {maximumAbsolute:0,maximumRelative:null});
});
test('scientifically different endpoints cannot produce timing ratios',()=>{
  const a={cases:{case:{...measurement,raw_endpoints:{density:[10,20]}}}};
  const b=structuredClone(a);
  assert.equal(report.compatible(a,b,'case'),true);
  b.cases.case.raw_endpoints.density[1]=21;
  assert.equal(report.compatible(a,b,'case'),false);
  assert.match(report.comparisonReason(a,b,'case'),/Different endpoint density/);
});
