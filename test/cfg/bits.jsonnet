local anodes = import "multi/anodes.jsonnet";
{
    xinsource : {
        type: "XinFileSource",
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
    
    xinsink: {
        type: "XinFileSink",
        data: {
            rebin: 6,
            filename: std.extVar("output"),
        },
    }

}
