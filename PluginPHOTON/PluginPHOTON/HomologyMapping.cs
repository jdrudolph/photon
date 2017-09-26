using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Text;
using BaseLib.Graphic;
using BaseLibS.Graph;
using BaseLibS.Param;
using PerseusApi.Matrix;
using Path = System.IO.Path;

namespace PluginPHOTON
{
    public class HomologyMapping : PluginInterop.Python.MatrixProcessing
    {
        public override string Heading => "Homology";
        public override string Name => "Map between humaan and mouse genes";
        public override string Description => "Map between mouse entrez gene ids and human entrez gene ids by homology from MGI.";
        public override bool HasButton => false;

        protected override string[] ReqiredPythonPackages => new[] { "perseuspy" };


        protected override string GetCodeFile(Parameters param)
        {
            byte[] code = (byte[])Properties.Resources.ResourceManager.GetObject("homology_mapping");
            var codeFile = Path.GetTempFileName();
            File.WriteAllText(codeFile, Encoding.UTF8.GetString(code));
            return codeFile;
        }


        protected override Parameter[] SpecificParameters(IMatrixData mdata, ref string errString)
        {
            if (mdata.StringColumnCount == 0)
            {
                errString = "Matrix does not have any text columns. Please add a 'GeneID' text column";
                return null;
            }
            var geneid = mdata.StringColumnNames.FindIndex(col => col.ToLower().Equals("geneid"));
            return new Parameter[] {
            new SingleChoiceParam("GeneID", Math.Max(0, geneid))
            {
                Values = mdata.StringColumnNames
            },
            new SingleChoiceParam("Map from", 0)
            {
                Values = new List<string> {"mouse to human", "human to mouse"}
            },
            };
        }
        
    }
}