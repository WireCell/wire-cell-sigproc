local anodes = import "multi/anodes.jsonnet";
{
    magnifysource : {
        type: "MagnifySource",
        data: {
            filename: std.extVar("input"),
            histtype: "raw",
        }
    },

    fieldresponse : {
        type: "FieldResponse",
        data: {
            filename: anodes.nominal.data.fields,
        }
    },

    perchanresp : {
        type : "PerChannelResponse",
        data : {
            filename: "calib_resp_v1.json.bz2",
        },
    },
    
    magnifysink: {
        type: "MagnifySink",
        data: {
            rebin: 6,
            // fixme: giving an input file here is evil.
            input_filename: std.extVar("input"),
            output_filename: std.extVar("output"),
        },
    }

}
