Instructor:Find more at https://karpathy.ai/ and https://x.com/karpathy

1,512,511 views  Feb 5, 2025
#This is a general audience deep dive into the Large Language Model (LLM) AI technology that powers ChatGPT and 
#related products. It is covers the full training stack of how the models are developed, along with mental models 
#of how to think about their "psychology", and how to get the best use them in practical applications.
# I have one "Intro to LLMs" video already from ~year ago, but that is just a re-recording of a random talk, 
# so I wanted to loop around and do a lot more comprehensive version.


#Instructor
#Andrej was a founding member at OpenAI (2015) and then Sr. Director of AI at Tesla (2017-2022), 
# and is now a founder at Eureka Labs, which is building an AI-native school. 
# His goal in this video is to raise knowledge and understanding of the state of the art in AI, 
# and empower people to effectively use the latest and greatest in their work.
Find more at https://karpathy.ai/ and https://x.com/karpathy
githib: https://github.com/karpathy

#Chatgpt

chaggpt_summary by the same author: https://youtu.be/7xTGNNLPyMI (Deep Dive into LLMs like ChatGPT)
#This is a general audience deep dive into the Large Language Model (LLM) AI technology that powers ChatGPT and related products. 
#It is covers the full training stack of how the models are developed, along with mental models of how to think about their "psychology", 
# and how to get the best use them in practical applications. I have one "Intro to LLMs" video already from ~year ago, but that is just 
# a re-recording of a random talk, so I wanted to loop around and do a lot more comprehensive version.


Chapters
00:00:00 introduction
00:01:00 pretraining data (internet)
00:07:47 tokenization
00:14:27 neural network I/O
00:20:11 neural network internals
00:26:01 inference
00:31:09 GPT-2: training and inference
00:42:52 Llama 3.1 base model inference
00:59:23 pretraining to post-training
01:01:06 post-training data (conversations)
01:20:32 hallucinations, tool use, knowledge/working memory
01:41:46 knowledge of self
01:46:56 models need tokens to think
02:01:11 tokenization revisited: models struggle with spelling
02:04:53 jagged intelligence
02:07:28 supervised finetuning to reinforcement learning
02:14:42 reinforcement learning
02:27:47 DeepSeek-R1
02:42:07 AlphaGo
02:48:26 reinforcement learning from human feedback (RLHF)
03:09:39 preview of things to come
03:15:15 keeping track of LLMs
03:18:34 where to find LLMs
03:21:46 grand summary

YouTube video link: https://youtu.be/EWvNQjAaOHw
#New 2h11m YouTube video: How I Use LLMs

#This video continues my general audience series. 
# The last one focused on how LLMs are trained, 
# so I wanted to follow up with a more practical guide of the entire LLM ecosystem,
# including lots of examples of use in my own life.

#Chapters give a sense of content:
##00:02:54 ChatGPT interaction under the hood
#00:13:12 Basic LLM interactions examples*
#00:18:03 Be aware of the model you're using, pricing tiers
#00:22:54 Thinking models and when to use them
#00:31:00 Tool use: internet search
#00:42:04 Tool use: deep research
#00:50:57 File uploads, adding documents to context
#00:59:00 Tool use: python interpreter, messiness of the ecosystem
#01:04:35 ChatGPT Advanced Data Analysis, figures, plots
#01:09:00 Claude Artifacts, apps, diagrams
#01:14:02 Cursor: Composer, writing code
#01:22:28 Audio (Speech) Input/Output
#01:27:37 Advanced Voice Mode aka true audio inside the model
#01:37:09 NotebookLM, podcast generation
#01:40:20 Image input, OCR
#01:47:02 Image output, DALL-E, Ideogram, etc.
#01:49:14 Video input, point and talk on app
#01:52:23 Video output, Sora, Veo 2, etc etc.
#01:53:29 ChatGPT memory, custom instructions
#01:58:38 Custom GPTs
#02:06:30 Summary



# + Excalidraw board we built up as notes also here as an image for an overview (and download link in the video description)
https://youtu.be/KJtZARuO3JY (Visualizing transformers and attention | Talk for TNG Big Tech Day '24; 
# thanks to Grant Sanderson)
#huggingface provides the detssils of the different models and the source code: https://huggingface.co/models
Fineweb and CommonCrawl from Hugging face is crucial in keeping data for training Models:https://huggingface.co/spaces/HuggingFaceFW/blogpost-fineweb-v1
tokenizer app :https://tiktokenizer.vercel.app/
excalidraw.com :https://excalidraw.com/
vizualize neeral network: https://bbycroft.net/llm
#Rent computation: Lambda Compute: https://aws.amazon.com/lambda/
https://lambdalabs.com/talk-to-an-engineer?primary_product_interest=Blackwell
https://cloud.lambdalabs.com/instances

#Other materials:
Links
ChatGPT https://chatgpt.com/
FineWeb (pretraining dataset): https://huggingface.co/spaces/Hugging...
Tiktokenizer: https://tiktokenizer.vercel.app/
Transformer Neural Net 3D visualizer: https://bbycroft.net/llm
llm.c Let's Reproduce GPT-2 https://github.com/karpathy/llm.c/dis...
Llama 3 paper from Meta: https://arxiv.org/abs/2407.21783
Hyperbolic, for inference of base model: https://app.hyperbolic.xyz/
InstructGPT paper on SFT: https://arxiv.org/abs/2203.02155
HuggingFace inference playground: https://huggingface.co/spaces/hugging...
DeepSeek-R1 paper: https://arxiv.org/abs/2501.12948
TogetherAI Playground for open model inference: https://api.together.xyz/playground
AlphaGo paper (PDF): https://discovery.ucl.ac.uk/id/eprint...
AlphaGo Move 37 video:    • Lee Sedol vs AlphaGo  Move 37 reactio...  
LM Arena for model rankings: https://lmarena.ai/
AI News Newsletter: https://buttondown.com/ainews
LMStudio for local inference https://lmstudio.ai/

The visualization UI I was using in the video: https://excalidraw.com/
The specific file of Excalidraw we built up: https://drive.google.com/file/d/1EZh5...
Discord channel for Eureka Labs and this video:   / discord  


openai gpt-2 model: 
https://github.com/openai/gpt-2/tree/master

Hyperbolic API: website for testin the models
https://docs.hyperbolic.xyz/docs/getting-started




