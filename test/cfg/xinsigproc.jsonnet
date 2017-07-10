// Example jsonnet command line: 
// $ jsonnet -J cfg -V detector=uboone -V input=foo.root sio/test/cfg/xinsigproc.jsonnet
// Similar wire-cell command line.

local params = import "params/chooser.jsonnet";
local wc = import "wirecell.jsonnet";
local anodes = import "multi/anodes.jsonnet";
local bits = import "bits.jsonnet";
local filters = import "filters.jsonnet";
local omni = import "omni.jsonnet";
[

    bits.magnifysource,

    anodes.nominal,

    bits.fieldresponse,

] + filters + [

    bits.perchanresp,

    bits.magnifysink,

    // omni.noisefilter,
    omni.sigproc,

    {
        type: "Omnibus",
        data: {
            source: wc.tn(bits.magnifysource),
            sink: wc.tn(bits.magnifysink),
            filters: [wc.tn(omni.sigproc)],
        }
    },


]
